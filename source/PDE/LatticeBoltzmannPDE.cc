// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <stddef.h>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <math.h>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>


#include "CoupledPDE/pdecoupling.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/Logging/log.hpp"
#include "DataInOut/programOptions.hh"
#include "DataInOut/resultHandler.hh"
#include "Domain/ansatzFct.hh"
#include "Domain/elem.hh"
#include "Domain/domain.hh"
#include "Domain/entityList.hh"
#include "Domain/grid.hh"
#include "Domain/resultInfo.hh"
#include "Driver/assemble.hh"
#include "Driver/basedriver.hh"
#include "Driver/formsContext.hh"
#include "Driver/solveStepSmooth.hh"
#include "Elements/basefe.hh"
#include "Forms/baseForm.hh"
#include "Forms/LatticeBoltzmannInt.hh"
#include "General/exception.hh"
#include "PDE/SinglePDE.hh"
#include "PDE/StdPDE.hh"
#include "PDE/eqnMap.hh"
#include "PDE/timestepping.hh"
#include "PDE/pseudoTS.hh"
#include "PDE/NonFEM/LatticeBoltzmann.hh"
#include "Utils/StdVector.hh"
#include "Utils/baseelemstoresol.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/TransferFunction.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Utils/result.hh"
#include "MatVec/vector.hh"
#include "MatVec/stdmatrix.hh"
#include "MatVec/crs_matrix.hh"
#include "OLAS/algsys/basesystem.hh"
#include "MatVec/stdmatrix.hh"
#include "LatticeBoltzmannPDE.hh"


namespace CoupledField {

using boost::numeric::ublas::matrix;


class BaseMaterial;
class SingleVector;

DECLARE_LOG(lbm_pde)
DEFINE_LOG(lbm_pde, "lbm_pde")

void test_matrix()
{
  compressed_matrix<double> cm(6,6);
  cm(0,0) = 10.0;
  cm(0,4) = -2;
  cm(1,0) = 3;
  cm(1,1) = 9;
  cm(1,5) = 3;
  cm(2,1) = 7;
  cm(2,2) = 8;
  cm(2,3) = 7;
  cm(3,0) = 3;
  cm(3,2) = 8;
  cm(3,3) = 7;
  cm(3,4) = 5;
  cm(4,1) = 8;
  cm(4,3) = 9;
  cm(4,4) = 9;
  cm(4,5) = 13;
  cm(5,1) = 4;
  cm(5,4) = 2;
  cm(5,5) = -1;

  compressed_matrix<double>::const_iterator1 iter;
//  iter = cm.find1(0, 1, 0);
//  std::cout << *iter << std::endl;
//  iter = cm.find1(0, 2, 0);
//  std::cout << *iter << std::endl;
//  iter = cm.find1(1, 1, 0);
//  std::cout << *iter << std::endl;
//  iter = cm.find1(1, 2, 0);
//  std::cout << *iter << std::endl;

  for(compressed_matrix<double>::const_iterator1 it1 = cm.begin1(); it1 != cm.end1(); ++it1) {
    for(compressed_matrix<double>::const_iterator2 it2 = it1.begin(); it2 != it1.end(); ++it2) {
      std::cout << cm(it2.index1(),it2.index2()) << std::endl;
    }
    std::cout << std::endl;
  }
//  std::cout<< cm << std::endl;
//  compressed_matrix<double, boost::numeric::ublas::column_major> ccs;
//  ccs = trans(cm);
//  std::cout << ccs << std::endl;

//  for(unsigned int i=0; i< cm.filled2();i++)
//    std::cout << "cm filled2: i=" << i << " -> index2_data=" << cm.index2_data()[i] << " v=" << cm.value_data()[i] << "\n";
//
//  // Create row array
//  for(unsigned int i=0;i< cm.filled1();i++)
//    std::cout << "cm filled1: i=" << i << " -> index1_data=" << cm.index1_data()[i] << " v=" << cm.value_data()[i] << "\n";
//
//  // this is the iterator over the rows
//  for(compressed_matrix<double>::const_iterator1 it = cm.begin1(); it != cm.end1(); ++it)
//  {
//    // it.i2 == 0
//    std::cout << "it1.i1=" << it.index1() << " it.i2=" << it.index2() << " it=" << *it << "\n";
//    // this is the iterator over the columns of the current row
//    for(compressed_matrix<double>::const_iterator2 it2 = it.begin(); it2 != it.end(); ++it2)
//       std::cout << "it2.i1=" << it2.index1() << " it2.i2=" << it2.index2() << " it2=" << *it2 << " == " <<  cm(it2.index1(), it2.index2()) << "\n" ;
//  }

  EXCEPTION("test_matrix");
}


LatticeBoltzmannPDE::LatticeBoltzmannPDE(Grid* grid, PtrParamNode pn) : SinglePDE(grid, pn)
{
  pdename_ = "LatticeBoltzmann";
  pdematerialclass_ = MECHANIC;
  firstTurn_ = true;
  nonLin_ = false;
  lbm = NULL;
  numWriteResults_ = 0;

  method_ = "mechanic";
  dirty_ = true; // not solved yet!

//  test_matrix();

  // LBM parametes
  omega_       = myParam_->Get("LBM/omega")->As<double>();
  maxWallTime_ = myParam_->Get("LBM/maxWallTime")->As<double>();
  maxIter_     = myParam_->Get("LBM/maxIter")->As<unsigned int>();
  convergence_ = myParam_->Get("LBM/convergence")->As<double>();
  writeFrequency_ = myParam_->Get("LBM/writeFrequency")->As<unsigned int>();
  bool plot    = myParam_->Get("LBM/plot")->As<bool>();

  PtrParamNode bcsl = myParam_->Get("bcsAndLoads");
  // only one inlet region
  u_x_ = bcsl->HasByVal("inlet", "dof", "x") ? bcsl->GetByVal("inlet", "dof", "x")->Get("value")->As<double>() : 0.0;
  u_y_ = bcsl->HasByVal("inlet", "dof", "y") ? bcsl->GetByVal("inlet", "dof", "y")->Get("value")->As<double>() : 0.0;
  u_z_ = bcsl->HasByVal("inlet", "dof", "z") ? bcsl->GetByVal("inlet", "dof", "z")->Get("value")->As<double>() : 0.0;

  iface.SetName("LatticeBoltzmannPDE::Iface");
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

  // n_q gives number of discrete velocities of LBM model
  if (domain->GetGrid()->GetDim() == 2)
    n_q_ = 9;
  else if (domain->GetGrid()->GetDim() == 3)
    n_q_ = 19;

  SetupElementMapping(grid);

  StdVector<Elem*> tmp;
  grid->GetElemsByName(tmp,"inlet");
  inlet.Resize(tmp.GetSize());
  for(unsigned int i = 0; i < inlet.GetSize(); ++i)
    inlet[i] = elem_to_idx[tmp[i]->elemNum];

  grid->GetElemsByName(tmp,"outlet");
  outlet.Resize(tmp.GetSize());
  for(unsigned int i = 0; i < outlet.GetSize(); ++i)
    outlet[i] = elem_to_idx[tmp[i]->elemNum];

  //Initializing storage for PDFs
  pdfs.Resize(n_elems * n_q_);

  if(iface_ == INTERNAL) {
    lbm = new LatticeBoltzmann(dim_, n_x_, n_y_, n_z_, u_x_, u_y_, u_z_, omega_, maxIter_, convergence_, plot, writeFrequency_);
  }

  velocityDirections = lbm->GetVelocityDirections();
  invDirections = lbm->GetInverseDirections();
}

LatticeBoltzmannPDE::~LatticeBoltzmannPDE()
{
  if(lbm) { delete lbm; lbm = NULL; }

}

void LatticeBoltzmannPDE::InitRegions(PtrParamNode pn, Grid* grid)
{
  // find the unique boundary region id
  ParamNodeList regions = pn->Get("regionList")->GetChildren();
  if(regions.GetSize() < 2)
    throw Exception("LatticeBoltzmann requires at least two regions where one has the boundary attribute set");

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
    throw Exception("LatticeBoltzmann requires a region with boundary attribute set");
}

void LatticeBoltzmannPDE::SetupElementMapping(Grid* grid)
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

void LatticeBoltzmannPDE::DefineIntegrators()
{
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
    BaseForm * bilinearStiff = new LatticeBoltzmannInt(actSDMat);

    BiLinFormContext * actIntDescrStiff = new BiLinFormContext(bilinearStiff, STIFFNESS );

    actIntDescrStiff->SetPtPdes(this, this);
    actIntDescrStiff->SetResults( results_[0], results_[0], actSDList, actSDList );
    assemble_->AddBiLinearForm( actIntDescrStiff );
    // Give entities and result to equation numbering class
    // and solution class
    eqnMap_->AddResult( *results_[0], actSDList );
  }
}

void LatticeBoltzmannPDE::InitCoupling(PDECoupling * coupling)
{
}

void LatticeBoltzmannPDE::DefineSolveStep()
{
  solveStep_ = new StdSolveStep(*this);
}
// returns number of iterations needed in LBM calculation until convergence
int LatticeBoltzmannPDE::GetNumWriteResults()
{
  return numWriteResults_;
}

void LatticeBoltzmannPDE::Solve()
{
  SetDirty(false); // assume the calculation will succeed. By this we may call WriteResults() which will call CalcResuls() which will not recursively call Salve() when dirty is set to false.

//  ResultHandler* rh = NULL;

//  if (writeFrequency_ > 1) {
//    rh = domain->GetResultHandler();
//    unsigned int mss = domain->GetDriver()->GetActSequenceStep();
//    // max steps is high. The number is only relevant for hdf5, but there a hard limit
//    rh->BeginMultiSequenceStep(mss, BasePDE::TRANSIENT, 9999);
//  }

  // infoNode_ is not set yet in the constructor
  PtrParamNode in = infoNode_->Get(ParamNode::HEADER)->Get("LBM");
  in->Get("omega")->SetValue(omega_);
  in->Get("maxIter")->SetValue(maxIter_);
  in->Get("maxWallTime")->SetValue(maxWallTime_);
  in->Get("convergence")->SetValue(convergence_);
  in->Get("iface")->SetValue(iface.ToString(iface_));

  // in the constructor we don't have the densities yet
  SetupElements();

  if(non_sing.IsEmpty()) // only once but we need elements
    SetNonSingualrityIndices();

  Timer timer;
  in = infoNode_->Get(ParamNode::PROCESS)->Get("stateProblem", progOpts->DoDetailedInfo() ? ParamNode::APPEND : ParamNode::INSERT);
  timer.Start(); // local timer
  state_.Start(); // global timer

  //for debugging
//  ExportMultipleFiles(elements);


  switch(iface_)
  {
  case EXT_MATLAB:
  case EXT_CFSxLBM:
  {

    executable = myParam_->Get("LBM")->Get("lbm")->As<std::string>();

    if(iface_ == EXT_CFSxLBM)
      ExportCFS2LBM(elements);
    else
      ExportMultipleFiles(elements);

    if(!boost::filesystem::exists(executable))
      EXCEPTION("Could not find executable '" + executable + "', might be not in path");

    std::cout << "++ Calling external LBM solver .. \n" << std::endl;
    std::string command = executable + " -o LBM2CFS.dat CFS2LBM.dat";
    int err = system(command.c_str());
    if (err)
      EXCEPTION("LBM simulation failed, no outputs available! \n");

    if(iface_ == EXT_CFSxLBM)
      ReadProbabilityDistribution("LBM2CFS.dat");
    else
      ReadProbabilityDistribution("node_steady.dat");

    break;
  }
  case INTERNAL:
  {
    LOG_DBG(lbm_pde) << "call internal LBM";
    LOG_DBG2(lbm_pde) << "elements\n" << ToString(elements, true, true);

    // this data is the final simulation but it lacks an additional propagation step for the adjoint calculation
    // We need to perform an additional propagation step as base for the adjoint system setup
    StdVector<double>* tmp = lbm->Iterate(elements, in->Get("LBM"));

    pdfs = *tmp;
    in->Get("original_pressure_drop")->SetValue(CalcPressureDrop());

    pdfs = *tmp;
    lbm->prop_step();

    pdfs = *tmp;
    in->Get("prop_step_pressure_drop")->SetValue(CalcPressureDrop());

    numWriteResults_ = lbm->GetNumWriteResults();

    break;
  }
  }

  timer.Stop();
  state_.Stop();
  in->Get("timer/cpu")->SetValue(timer.GetCPUTime());
  in->Get("timer/wall")->SetValue(timer.GetWallTime());

  in = infoNode_->Get(ParamNode::SUMMARY)->Get("stateProblem");
  in->Get("totalTimer/cpu")->SetValue(state_.GetCPUTime());
  in->Get("totalTimer/wall")->SetValue(state_.GetWallTime());
  in->Get("totalTimer/calls")->SetValue(state_.GetCalls());

  // dirty_ = false; already set to false in the beginning
}

void LatticeBoltzmannPDE::SetupSensitivityAnalysis(StdVector<double>& ux, StdVector<double>& uy, StdVector<double>& dloc, StdVector<double>& weights)
{
  ux.Resize(n_elems);
  uy.Resize(n_elems);
  dloc.Resize(n_elems);

  // macroscopic values for each m_node
  for(unsigned int i = 0; i < n_elems; i++)
  {
    double density = CalcLBMDensity(i);
    ux[i] = CalcVelocityX(i, density);
    uy[i] = CalcVelocityY(i, density);
    dloc[i] = density;

    // sets quadrature weights
    weights.Resize(n_q_);
    weights[0] = 4./9.;
    for (int i=1; i<5; i++)
      weights[i] = 1./9.;
    for (int i=5; i<9; i++)
      weights[i] = 1./36.;
    }
  }

}


void LatticeBoltzmannPDE::SensitivityAnalysis(TransferFunction* tf, Function* f, DesignSpace* space)
{
  // needs to be called after Solve() !!

  Timer timer;
  timer.Start();
  adjoint_.Start();

  // Initialization of the data structures
  StdVector<double> ux(n_elems);
  StdVector<double> uy(n_elems);
  StdVector<double> dloc(n_elems);
  StdVector<double> weights(n_q_);

  SetupSensitivityAnalysis(ux, uy, dloc, weights);

  // derivative of the residual with respect to design variables
  Vector<double> dRds(n_q_ * n_elems);
  dRds.Init(0.0);

  // derivative of the collision operator with respect to p
  mapped_matrix<double> col_jacobi(n_q_ * n_elems, n_q_ * n_elems);
  //compressed_matrix<double> col_jacobi(9 * n_elems, 9 * n_elems);

  StdVector<double> dfeqdux(n_q_);
  StdVector<double> dfeqduy(n_q_);
  Matrix<double> block(n_q_,n_q_);
  //cout<<"sens start correct"<<endl;
  for(unsigned int index = 0; index < n_elems ; index++)
  {
    block.InitValue(0.0);
    // interior
    if(elements[index] >= 0)  {
      // Jacobi matrix of the collision operator with respect to the design variables
      d_collision_step_d_f(index, block, dfeqdux, dfeqduy, ux, uy, dloc, weights);

      // simple choice of the mapping between p and s
      DesignElement* de = space->Find(idx_to_elem[index], DesignElement::DENSITY);
      // scale = -penal * pow(por[0][index],(penal-1.));
      double scale = -1.0 * tf->Derivative(de, DesignElement::SMART);
      for (unsigned int k=0;k<n_q_;k++)
        dRds[index * n_q_ +k] = omega_ * (dfeqdux[k] * scale * ux[index] + dfeqduy[k] * scale * uy[index]);
    }
    else if(elements[index] == LBM_NODE_TYPE_BB) {
      d_bounceback_d_f(block); // bounce-back sensitivities
    }
    else if(elements[index] == LBM_NODE_TYPE_INLET) {
      d_inflow_d_f(index,block, weights); // inlet derivative with respect to f
    }
    else if(elements[index] == LBM_NODE_TYPE_OUTLET) {
      d_outflow_d_f(index,block, ux, uy, dloc, weights); // outlet derivative with respect to f
    }
    else {
      assert(false);
    }
    // fill transpose of block in col_jacobi
    for(unsigned int k = 0; k <n_q_; k++)
      for(unsigned int l = 0; l<n_q_; l++)
        col_jacobi(index * n_q_ + k,index * n_q_ + l) = block[l][k];
  }
  
  double d_collision_setup = timer.GetCPUTime();

  // the real jacobi combines d_collision with d_propagation
  mapped_matrix<double> Jacobi(n_q_ * n_elems , n_q_ * n_elems);
  // compressed_matrix<double> Jacobi(9 * n_elems , 9 * n_elems);
  d_propagate_d_f(Jacobi,col_jacobi);

  for(unsigned int i=0;i < n_q_ * n_elems; i++)
    Jacobi(i,i) -= 1.;

  double d_propagation_setup = timer.GetCPUTime(); 

  //Data for Pardiso solver
  compressed_matrix<double> Jacobi_new(n_elems * n_q_, n_elems * n_q_, Jacobi.nnz());

  // Delete singular rows from the Jacobian
  DeleteSingularities(Jacobi,Jacobi_new);

  double delete_sing_setup = timer.GetCPUTime();

  // right-hand side
  Vector<double> b = d_pressuredrop_d_f(ux, uy, dloc);

  double rhs_setup = timer.GetCPUTime();
   
  // use the cfssystem matrix to solve the adjoint system
  StdMatrix* mat = algsys_->GetSysMat();
  LOG_DBG(lbm_pde) << "SA: mat structure=" << mat->GetStructureType();
  LOG_DBG(lbm_pde) << "SA: mat storage=" << mat->GetStorageType();

  CRS_Matrix<double>* crs = dynamic_cast<CRS_Matrix<double>*>(mat);
  assert(crs != NULL);

  crs->SetSize(n_elems * n_q_, n_elems * n_q_, Jacobi_new.nnz());
  matrix_sparse_to_crs(Jacobi_new, crs->GetDataPointer(), crs->GetRowPointer(), crs->GetColPointer());

  // time to setup adjoint system before solving
  double setup_wall = timer.GetWallTime();
  double setup_cpu = timer.GetCPUTime();
  
  algsys_->InitRHS(b);

  PtrParamNode analysis_id = domain->GetDriver()->CreateAnalysisId("lbm_adjoint", 0);
  algsys_->SetupSolver(analysis_id);
  // !!!!! The problem needs to be solved transposed!!! For pardiso this is <Transposed>transposed</Transposed>!!!
  algsys_->Solve(analysis_id);
  Vector<double> sol;
  algsys_->GetSolutionVal(sol);

  for(unsigned int e = 0; e < f->elements.GetSize(); e++)
  {
    DesignElement* de = f->elements[e];
    unsigned int idx = elem_to_idx[de->elem->elemNum]; // lbm idx
    double val = -1.0 * sol.Inner(dRds, idx * n_q_, (idx + 1) * n_q_);
    de->AddGradient(f, val);
  }

  timer.Stop();
  adjoint_.Stop();
  PtrParamNode adjoint = infoNode_->Get(ParamNode::PROCESS)->Get("adjoint", progOpts->DoDetailedInfo() ? ParamNode::APPEND : ParamNode::INSERT);
  adjoint->Get("timer/cpu")->SetValue(timer.GetCPUTime());
  adjoint->Get("timer/wall")->SetValue(timer.GetWallTime());
  adjoint->Get("setupTimer/cpu")->SetValue(setup_cpu);
  adjoint->Get("setupTimer/wall")->SetValue(setup_wall);
  adjoint->Get("setupTimer/d_coll")->SetValue(d_collision_setup);
  adjoint->Get("setupTimer/d_prop")->SetValue(d_propagation_setup - d_collision_setup);
  adjoint->Get("setupTimer/del_sing")->SetValue(rhs_setup - delete_sing_setup);
  adjoint->Get("setupTimer/rhs")->SetValue(delete_sing_setup - d_propagation_setup);
  

  adjoint = infoNode_->Get(ParamNode::SUMMARY)->Get("adjoint");
  adjoint->Get("totalTimer/cpu")->SetValue(adjoint_.GetCPUTime());
  adjoint->Get("totalTimer/wall")->SetValue(adjoint_.GetWallTime());
  adjoint->Get("totalTimer/calls")->SetValue(adjoint_.GetCalls());
}

void LatticeBoltzmannPDE::matrix_sparse_to_crs(compressed_matrix<double>& M, double* a, unsigned int* ia, unsigned int* ja)
{
  //conversion of sparse matrix to compressed row storage format
  // Create value array and column array
  for(unsigned int i=0;i<M.filled2();i++) {
    a[i] = M.value_data()[i];
    ja[i] = M.index2_data()[i];
  }
  // Create row array
  for(unsigned int i=0;i<M.filled1();i++) {
    ia[i] = M.index1_data()[i];
  }
}


//void LatticeBoltzmannPDE::matrix_sparse_to_crs(compressed_matrix<double>& M, double* a, unsigned int* ia, unsigned int* ja)
void LatticeBoltzmannPDE::matrix_sparse_to_crs(mapped_matrix<double>& M, double* a, unsigned int* ia, unsigned int* ja)
{
  /*
  //conversion of sparse matrix to compressed row storage format
  // Create value array and column array
  for(unsigned int i=0;i<M.filled2();i++) {
    a[i] = M.value_data()[i];
    ja[i] = M.index2_data()[i];
  }
  // Create row array
  for(unsigned int i=0;i<M.filled1();i++) {
    ia[i] = M.index1_data()[i];
  }
  */
}



// void LatticeBoltzmannPDE::DeleteSingularities(const compressed_matrix<double> & M,compressed_matrix<double> & output)
void LatticeBoltzmannPDE::DeleteSingularities(const mapped_matrix<double> & M, compressed_matrix<double> & output)
{
  std::set<unsigned int> sing;
  compressed_matrix<double> test(n_elems * n_q_, n_elems * n_q_, output.nnz());

  for(unsigned int i = 0, n = non_sing.GetSize(); i < n; i++)
    sing.insert(non_sing[i] + 1);

  // delete singular rows of the Jacobian in the sensitivity analysis for the optimization
  double val;
  std::set<unsigned int>::iterator sing_it1, sing_it2;
  // for(compressed_matrix<double>::const_iterator1 it = M.begin1(); it != M.end1(); ++it) {
  for(mapped_matrix<double>::const_iterator1 it = M.begin1(); it != M.end1(); ++it) {
    // for(compressed_matrix<double>::const_iterator2 it2 = it.begin(); it2 != it.end(); ++it2) {
    for(mapped_matrix<double>::const_iterator2 it2 = it.begin(); it2 != it.end(); ++it2) {
      val = M(it2.index1(),it2.index2());
      if (abs(val) > 1e-20) {
        //it2.index1():row, it2.index2():column
        sing_it1 = sing.find(it2.index1()+1);
        sing_it2 = sing.find(it2.index2()+1);
        if (sing_it1 != sing.end() && sing_it2 != sing.end()) {
          // if elements found in sing
          //output(it2.index2(),it2.index1()) = val;
          output(it2.index1(),it2.index2()) = val;
          LOG_DBG3(lbm_pde) << "DS: o("<< it2.index2() << ", " << it2.index1() << ") = " << val;
        }
      }
    }
  }
  std::set<unsigned int>::iterator it = sing.begin();
  for(unsigned int i=0;i<output.size1();i++) {
    if (it == sing.end()) {
      output(i,i) = 1.;
      LOG_DBG3(lbm_pde) << "DS: a i=" << i;
    } else {
      if (i+1 != *it) {
        output(i,i) = 1.;
        LOG_DBG3(lbm_pde) << "DS: b i=" << i;
      } else {
        it++;
      }
    }
  }
}



void LatticeBoltzmannPDE::d_collision_step_d_f(unsigned int index, Matrix<double>& block, StdVector<double>& dfeqdux, StdVector<double>& dfeqduy, const StdVector<double>& ux, const StdVector<double>& uy, const StdVector<double>& dloc, const StdVector<double>& weight)
{
  // gradient of the collision step with respect to the design variables
  // partial derivative of f^eq with respect to rho, ux and uy: CORRECT (FDM)
  StdVector<double> dfeqdrho;
  dfeqdrho.Resize(n_q_);
//  double E[dim_][n_q_] = {{0, 1, 0, -1, 0, 1, -1, -1, 1},{0, 0, 1, 0, -1, 1, 1, -1, -1}};
  double scale = 1. - elements[index]; // 1.-pow(por[0][index],penal);
  assert(scale >= 0);
  LOG_DBG3(lbm_pde) << "d_collision_step_d_f: index = " << index << " scale=" << scale;
  //LOG_DBG3(lbm_pde) << "d_collision_step_d_f: index = " << index << " pdf=" << StdVector<double>::ToString(9, &pdfs.GetPointer()[index * 9]);
  //partial derivative of f^eq with respect to rho, ux and uy: CORRECT (FDM)
  double us1,us2,dot,norm;
  double const1 = 9. / 2.;
  double const2 = 3. / 2.;
  LatticeBoltzmann::PropTransform* transform;

  //gradient of u_x with respect to f
  StdVector<double> duxdf;
  duxdf.Resize(n_q_);
  //gradient of u_y with respect to f
  StdVector<double> duydf;
  duydf.Resize(n_q_);
//  double tmp[n_q_] = {0.,1.,0.,-1.,0.,1.,-1.,-1.,1.};
//  double tmp2[n_q_] = {0.,0.,1.,0.,-1.,1.,1.,-1.,-1.};

  double invdloc = 1.0 / dloc[index];

  for (unsigned int i = 0; i < n_q_; i++)
  {
    transform = &(*velocityDirections)[i];
    us1 = scale * ux[index];
    us2 = scale * uy[index];
    dot = transform->off_x * us1 + transform->off_y * us2;
    norm = us1 * us1 + us2 * us2;
    dfeqdrho[i] = weight[i] * (1. + 3.*dot + const1*dot*dot - const2* norm);
    dfeqdux[i] = weight[i] * dloc[index]*(3. * transform->off_x + 9. * transform->off_x * dot - 3. * us1);
    dfeqduy[i] = weight[i] * dloc[index]*(3. * transform->off_y + 9. * transform->off_y * dot - 3. * us2);
    //LOG_DBG3(lbm_pde) << "d_collision_step_d_f: w = " << weight[i] << " dloc=" << dloc[index] << " dot=" << dot;

    duxdf[i] = scale*invdloc*(-ux[index] + transform->off_x);
    duydf[i] = scale*invdloc*(-uy[index] + transform->off_y);
  }

  // partial derivatives of f^eq with respect to f: CORRECT (FDM)
  StdVector< StdVector<double> > dfeqdf;
  dfeqdf.Resize(n_q_);
  for (unsigned int dir = 0; dir < n_q_; dir++ ){
    dfeqdf[dir].Resize(n_q_);
  }

  for (unsigned int i=0; i<n_q_; i++)
    for (unsigned int j= 0; j<n_q_; j++)
//      dfeqdf[i][j] = dfeqdrho[i] * drhodf + dfeqdux[i]*duxdf[j] + dfeqduy[i]*duydf[j]; NOTE: drhodf = [1,1,1,1,1,1,1,1,1]
      dfeqdf[i][j] = dfeqdrho[i] + dfeqdux[i]*duxdf[j] + dfeqduy[i]*duydf[j];

  //LOG_DBG3(lbm_pde) << "d_collision_step_d_f: index = " << index << " dfeqdux=" << dfeqdux.ToString() ;

  //partial derivative of collision operator with respect to f: CORRECT (FDM)
    //std::cout << "dcollision_step index=" << index << " block i=" << i << ": ";
  for (unsigned int i=0;i<n_q_;i++)
  {
    for (unsigned int j= 0;j<n_q_;j++)
    {

      //      if (i == j)
      //        block[i][j] = 1. - omega_ *(1.- dfeqdf[i][j]);
      //      else
      //        block[i][j] = omega_ * dfeqdf[i][j];
      //std::cout << block[i][j] << " ";
      block[i][j] = (i == j) * (1 - omega_) + omega_ * dfeqdf[i][j];
    }
    //std::cout << "\n";
  }

}

void LatticeBoltzmannPDE::d_bounceback_d_f(Matrix<double>& block)
{
  // bounce-back sensitivities; gradient is just a permutation matrix
  for (unsigned int dir = 0; dir < n_q_; dir++) {
    block[dir][(*invDirections)[dir]] = 1.0;
  }
}

void LatticeBoltzmannPDE::d_inflow_d_f(int index, Matrix<double>& block, StdVector<double>& weight)
{
  // gradient of the velocity inlet boundary with respect to the design variables
//  double dfeqdrho[9];
  StdVector<double> dfeqdrho;
  dfeqdrho.Resize(n_q_);
//  double E[2][9] = {{0, 1, 0, -1, 0, 1, -1, -1, 1},{0, 0, 1, 0, -1, 1, 1, -1, -1}};
  double us1 = 1.0 * u_x_;
  double us2 = 1.0 * u_y_;
  double dot, norm;
  LatticeBoltzmann::PropTransform* transform;
  for (unsigned int i=0; i < n_q_; i++)
  {
    transform = &(*velocityDirections)[i];
    dot = transform->off_x * us1 + transform->off_y * us2;
    norm = us1*us1+us2*us2;
    dfeqdrho[i] = weight[i]*(1. + 3*dot +4.5*dot*dot - 1.5* norm);
    // LOG_DBG3(lbm_pde) << "d_inflow_d_f index=" << index << " i=" << i << " us1=" << us1 << " us2=" << us2 << " dot=" << dot << " norm=" << norm;
  }
  // LOG_DBG3(lbm_pde) << "d_inflow_d_f index=" << index << " dfeqdrho=" << StdVector<double>::ToString(9, &dfeqdrho[0]);

  for (unsigned int i = 0; i < n_q_; i++)
    for (unsigned int j = 0; j < n_q_; j++)
      block[i][j] = dfeqdrho[i];

}

void LatticeBoltzmannPDE::d_outflow_d_f(int index, Matrix<double>& block, StdVector<double>& ux, StdVector<double>& uy, StdVector<double>& dloc, StdVector<double>& weight)
{
  // gradient of the density outlet boundary condition with respect to the design variables
//  double dfeqdux[9],dfeqduy[9];
  StdVector<double> dfeqdux, dfeqduy;
  dfeqdux.Resize(n_q_);
  dfeqduy.Resize(n_q_);
//  double E[2][9] = {{0, 1, 0, -1, 0, 1, -1, -1, 1},{0, 0, 1, 0, -1, 1, 1, -1, -1}};
  LatticeBoltzmann::PropTransform* transform;
  double scale = 1.; // no design
  double density = 1.0;  // we have no other case
  double us1,us2, dot;
  for (unsigned int i = 0; i < n_q_; i++)
  {
    transform = &(*velocityDirections)[i];
    us1 = scale*ux[index];
    us2 = scale*uy[index];
    dot = transform->off_x * us1 + transform->off_y * us2;
    dfeqdux[i] = weight[i]*density*(3 *  transform->off_x  + 9. * transform->off_x * dot - 3 * us1);
    dfeqduy[i] = weight[i]*density*(3 * transform->off_y + 9.* transform->off_y * dot - 3 * us2);
  }

//  double duxdf[9];
//  double duydf[9];
//  double tmp[9] = {0.,1.,0.,-1.,0.,1.,-1.,-1.,1.};
//  double tmp2[9] = {0.,0.,1.,0.,-1.,1.,1.,-1.,-1.};

  //gradient of u_x/u_y with respect to f
  StdVector<double> duxdf, duydf;
  duxdf.Resize(n_q_);
  duydf.Resize(n_q_);

  double invdloc = 1.0 / dloc[index];

  for (unsigned int i = 0; i < n_q_ ; i++) {
    transform = &(*velocityDirections)[i];
    duxdf[i] = invdloc * (-ux[index] + transform->off_x);
    duydf[i] = invdloc * (-uy[index] + transform->off_y);
  }

  //gradient of u_x with respect to f
  for (unsigned int i = 0; i < n_q_; i++)
    for (unsigned int j = 0;j < n_q_; j++)
      block[i][j] = dfeqdux[i]*duxdf[j] + dfeqduy[i]*duydf[j];

}

//void LatticeBoltzmannPDE::d_propagate_d_f_inDir(mapped_matrix<double>& Jprop, const mapped_matrix<double>& J, int dir) {
//  int rows1, rows2;
//  mapped_matrix<double>::const_iterator1 iter1, iter2;
//  switch (dir)
//  {
//    case 0:
//      rows1 = (y-1)*n_x_*9+(x-1)*9;
//      rows2 = 0;
//      iter1 = J.find1(0, rows1, 0);
//      iter2 = iter1.end();
//      break;
//    case 1:
//      rows1=((y-1)*(n_x_)*9+1);
//      rows2=((y-1)*(n_x_)*9+9+1);
//      iter1 = J.find1(0, rows1, 0);
//      iter2 = J.find1(0, rows2, 0);
//      break;
//  }
//  for(mapped_matrix<double>::const_iterator2 it = iter.begin(); it != iter.end(); ++it) {
//    Jprop(rows1,it.index2())=J(rows1,it.index2());
//  }
//  for(mapped_matrix<double>::const_iterator2 it = iter.begin(); it != iter.end(); ++it) {
//    Jprop(rows1,it.index2())=J(rows1,it.index2());
//  }
//}

// void LatticeBoltzmannPDE::d_propagate_d_f(compressed_matrix<double>& Jprop, const compressed_matrix<double>& J)
void LatticeBoltzmannPDE::d_propagate_d_f(mapped_matrix<double>& Jprop, const mapped_matrix<double>& J)
{
  //gradient of the propagations step with respect to the design variables by Georg Pingen, University of Colorado (CU), Boulder, Colorado
  // f0
  int rows1,rows2;
  // compressed_matrix<double>::const_iterator1 iter;
  mapped_matrix<double>::const_iterator1 iter;
  // fdist 0 */
  for(unsigned int y = 0;y < n_y_;y++) {
    for(unsigned int x = 0;x < n_x_;x++) {
      rows1 = GetPdfIndex(x,y,Q_0);
      iter = J.find1(0, rows1, 0);
      for(mapped_matrix<double>::const_iterator2 it = iter.begin(); it != iter.end(); ++it) {
        Jprop(rows1,it.index2())=J(rows1,it.index2());
        // std::cout << "dp f_0 (" << rows1 << ", " << it.index2() << ") <- (" << rows1 << ", " << it.index2() << ") : 1\n";
      }
    }
  }

  //* fdist 1 */
  for(unsigned int y = 0; y < n_y_; y++) {
    rows1 = GetPdfIndex(0,y,Q_E);
    rows2 = GetPdfIndex(1,y,Q_E);
    iter = J.find1(0, rows1, 0);
    for(mapped_matrix<double>::const_iterator2 it = iter.begin(); it != iter.end(); ++it) {
      Jprop(rows1,it.index2())= J(rows1,it.index2());
      // std::cout << "dp f_1 (" << rows1 << ", " << it.index2() << ") <- (" << rows1 << ", " << it.index2() << ") : 2\n";
    }
    iter = J.find1(0, rows2, 0);
    for(mapped_matrix<double>::const_iterator2 it = iter.begin(); it != iter.end(); ++it) {
      Jprop(rows1,it.index2()) = Jprop(rows1,it.index2()) + J(rows2,it.index2());
      // std::cout << "dp f_1 (" << rows1 << ", " << it.index2() << ") <- (" << rows2 << ", " << it.index2() << ") : 3\n";
    }

  }

  for(unsigned int y = 0; y < n_y_; y++) {
    for(unsigned int x = 1; x < (n_x_ - 1); x++) {
      rows1 = GetPdfIndex(x,y,Q_E);
      rows2 = GetPdfIndex(x+1,y,Q_E);
      iter = J.find1(0, rows2, 0);
      for(mapped_matrix<double>::const_iterator2 it = iter.begin(); it != iter.end(); ++it) {
        Jprop(rows1,it.index2())=J(rows2,it.index2());
        // std::cout << "dp f_1 (" << rows1 << ", " << it.index2() << ") <- (" << rows2 << ", " << it.index2() << ") : 4 \n";
      }
    }
  }

  //* fdist2 */
  for (unsigned int x = 0; x < n_x_; x++) {
    rows1 = GetPdfIndex(x,0,Q_N);
    rows2 = GetPdfIndex(x,1,Q_N);
    iter = J.find1(0, rows1, 0);
    for (mapped_matrix<double>::const_iterator2 it = iter.begin();
        it != iter.end(); ++it) {
      Jprop(rows1, it.index2()) = J(rows1, it.index2());
      // std::cout << "dp f_2 (" << rows1 << ", " << it.index2() << ") <- (" << rows1 << ", " << it.index2() << ") : 8 \n";
    }
    iter = J.find1(0, rows2, 0);
    for (mapped_matrix<double>::const_iterator2 it = iter.begin();
        it != iter.end(); ++it) {
      Jprop(rows1, it.index2()) = Jprop(rows1, it.index2())
                  + J(rows2, it.index2());
      // std::cout << "dp f_3 (" << rows1 << ", " << it.index2() << ") <- (" << rows2 << ", " << it.index2() << ") : 9 \n";
    }
  }

  for (unsigned int y = 1; y < n_y_ - 1; y++) {
    for (unsigned int x = 0; x < n_x_; x++) {
      rows1 = GetPdfIndex(x,y,Q_N);
      rows2 = GetPdfIndex(x,y+1,Q_N);
      iter = J.find1(0, rows2, 0);

      for (mapped_matrix<double>::const_iterator2 it = iter.begin();
          it != iter.end(); ++it) {
        Jprop(rows1, it.index2()) = J(rows2, it.index2());
        // std::cout << "dp f_3 (" << rows1 << ", " << it.index2() << ") <- (" << rows2 << ", " << it.index2() << ") : 10 \n";
      }
    }
  }

  //* fdist3 */
  for(unsigned int y = 0; y < n_y_; y++) {
    rows1 = GetPdfIndex(n_x_-1,y,Q_W);
    rows2 = GetPdfIndex(n_x_-2,y,Q_W);
    iter = J.find1(0, rows1, 0);
    for(mapped_matrix<double>::const_iterator2 it = iter.begin(); it != iter.end(); ++it) {
      Jprop(rows1,it.index2())= J(rows1,it.index2());
      // std::cout << "dp f_3 (" << rows1 << ", " << it.index2() << ") <- (" << rows1 << ", " << it.index2() << ") : 5 \n";
    }
    iter = J.find1(0, rows2, 0);
    for(mapped_matrix<double>::const_iterator2 it = iter.begin(); it != iter.end(); ++it) {
      Jprop(rows1,it.index2()) = Jprop(rows1,it.index2()) + J(rows2,it.index2());
      // std::cout << "dp f_3 (" << rows1 << ", " << it.index2() << ") <- (" << rows2 << ", " << it.index2() << ") : 6 \n";

    }
  }

  for(unsigned int y = 0; y < n_y_; y++) {
    for(unsigned int x = n_x_- 2; x >= 1;x--) {
      rows1 = GetPdfIndex(x,y,Q_W);
      rows2 = GetPdfIndex(x-1,y,Q_W);
      iter = J.find1(0, rows2, 0);
      for(mapped_matrix<double>::const_iterator2 it = iter.begin(); it != iter.end(); ++it) {
        Jprop(rows1,it.index2())=J(rows2,it.index2());
        // std::cout << "dp f_3 (" << rows1 << ", " << it.index2() << ") <- (" << rows2 << ", " << it.index2() << ") : 7 \n";
      }
    }
  }

  //* fdist4 */
  for(unsigned int x = 0; x < n_x_;x++) {
    rows1 = GetPdfIndex(x,n_y_-1,Q_S);
    rows2 = GetPdfIndex(x,n_y_-2,Q_S);
    iter = J.find1(0, rows1, 0);

    for(mapped_matrix<double>::const_iterator2 it = iter.begin(); it != iter.end(); ++it) {
      Jprop(rows1,it.index2())= J(rows1,it.index2());
      // std::cout << "dp f_3 (" << rows1 << ", " << it.index2() << ") <- (" << rows1 << ", " << it.index2() << ") : 11 \n";
    }
    iter = J.find1(0, rows2, 0);

    for(mapped_matrix<double>::const_iterator2 it = iter.begin(); it != iter.end(); ++it) {
      Jprop(rows1,it.index2()) = Jprop(rows1,it.index2()) + J(rows2,it.index2());
      // std::cout << "dp f_3 (" << rows1 << ", " << it.index2() << ") <- (" << rows2 << ", " << it.index2() << ") : 12 \n";
    }
  }

  for(unsigned int y = n_y_ - 2; y >= 1;y--) {
    for(unsigned int x = 0; x < n_x_;x++) {
      rows1 = GetPdfIndex(x,y,Q_S);
      rows2 = GetPdfIndex(x,y-1,Q_S);
      iter = J.find1(0, rows2, 0);

      for(mapped_matrix<double>::const_iterator2 it = iter.begin(); it != iter.end(); ++it) {
        Jprop(rows1,it.index2())=J(rows2,it.index2());
        // std::cout << "dp f_3 (" << rows1 << ", " << it.index2() << ") <- (" << rows2 << ", " << it.index2() << ") : 13 \n";
      }
    }
  }

  //* fdist5 */

  // left edge with lower left corner
  for(unsigned int y = 0;y < n_y_-1; y++) {
    rows1 = GetPdfIndex(0,y,Q_NE);
    rows2 = GetPdfIndex(1,y+1,Q_NE);
    iter = J.find1(0, rows1, 0);

    for(mapped_matrix<double>::const_iterator2 it = iter.begin(); it != iter.end(); ++it) {
      Jprop(rows1,it.index2())= J(rows1,it.index2());
      // std::cout << "dp f_5 (" << rows1 << ", " << it.index2() << ") <- (" << rows1 << ", " << it.index2() << ") : 16 \n";
    }
    iter = J.find1(0, rows2, 0);

    for(mapped_matrix<double>::const_iterator2 it = iter.begin(); it != iter.end(); ++it) {
      Jprop(rows1,it.index2()) = Jprop(rows1,it.index2()) + J(rows2,it.index2());
      // std::cout << "dp f_5 (" << rows1 << ", " << it.index2() << ") <- (" << rows2 << ", " << it.index2() << ") : 17 \n";
    }
  }

  // lower edge
  for(unsigned int x = 1; x < n_x_-1;x++) {
    rows1 = GetPdfIndex(x,0,Q_NE);
    rows2 = GetPdfIndex(x+1,1,Q_NE);
    iter = J.find1(0, rows1, 0);

    for(mapped_matrix<double>::const_iterator2 it = iter.begin(); it != iter.end(); ++it) {
      Jprop(rows1,it.index2())= J(rows1,it.index2());
      // std::cout << "dp f_5 (" << rows1 << ", " << it.index2() << ") <- (" << rows1 << ", " << it.index2() << ") : 18 \n";
    }
    iter = J.find1(0, rows2, 0);

    for(mapped_matrix<double>::const_iterator2 it = iter.begin(); it != iter.end(); ++it) {
      Jprop(rows1,it.index2()) = Jprop(rows1,it.index2()) + J(rows2,it.index2());
      // std::cout << "dp f_5 (" << rows1 << ", " << it.index2() << ") <- (" << rows2 << ", " << it.index2() << ") : 19 \n";
    }
  }

  // inner elements
  for(unsigned int y = 1; y < n_y_-1;y++) {
    for(unsigned int x = 1;x < n_x_-1;x++) {
      rows1 = GetPdfIndex(x,y,Q_NE);
      rows2 = GetPdfIndex(x+1,y+1,Q_NE);;
      iter = J.find1(0, rows2, 0);

      for(mapped_matrix<double>::const_iterator2 it = iter.begin(); it != iter.end(); ++it) {
        Jprop(rows1,it.index2())=J(rows2,it.index2());
        // std::cout << "dp f_5 (" << rows1 << ", " << it.index2() << ") <- (" << rows2 << ", " << it.index2() << ") : 20 \n";
      }
    }
  }


  //* fdist6 */

  // right edge with lower right corner
  for(unsigned int y = 0; y < n_y_-1; y++) {
    rows1 = GetPdfIndex(n_x_-1,y,Q_NW);
    rows2 = GetPdfIndex(n_x_-2,y+1,Q_NW);
    iter = J.find1(0, rows1, 0);

    for(mapped_matrix<double>::const_iterator2 it = iter.begin(); it != iter.end(); ++it) {
      Jprop(rows1,it.index2())= J(rows1,it.index2());
      // std::cout << "dp f_6 (" << rows1 << ", " << it.index2() << ") <- (" << rows1 << ", " << it.index2() << ") : 22 \n";

    }
    iter = J.find1(0, rows2, 0);

    for(mapped_matrix<double>::const_iterator2 it = iter.begin(); it != iter.end(); ++it) {
      Jprop(rows1,it.index2()) = Jprop(rows1,it.index2()) + J(rows2,it.index2());
      // std::cout << "dp f_6 (" << rows1 << ", " << it.index2() << ") <- (" << rows2 << ", " << it.index2() << ") : 23 \n";
    }
  }

  // lower edge without corners
  for(unsigned int x = n_x_ - 2; x >= 1; x--) {
    rows1 = GetPdfIndex(x,0,Q_NW);
    rows2 = GetPdfIndex(x-1,1,Q_NW);
    iter = J.find1(0, rows1, 0);

    for(mapped_matrix<double>::const_iterator2 it = iter.begin(); it != iter.end(); ++it) {
      Jprop(rows1,it.index2())= J(rows1,it.index2());
      // std::cout << "dp f_6 (" << rows1 << ", " << it.index2() << ") <- (" << rows1 << ", " << it.index2() << ") : 24 \n";
    }
    iter = J.find1(0, rows2, 0);

    for(mapped_matrix<double>::const_iterator2 it = iter.begin(); it != iter.end(); ++it) {
      Jprop(rows1,it.index2()) = Jprop(rows1,it.index2()) + J(rows2,it.index2());
      // std::cout << "dp f_6 (" << rows1 << ", " << it.index2() << ") <- (" << rows2 << ", " << it.index2() << ") : 25 \n";
    }
  }

  //inner elements
  for(unsigned int y = 1; y < n_y_-1; y++) {
    for(unsigned int x = n_x_ - 2; x >= 1; x--) {
      rows1 = GetPdfIndex(x,y,Q_NW);
      rows2 = GetPdfIndex(x-1,y+1,Q_NW);
      iter = J.find1(0, rows2, 0);

      for(mapped_matrix<double>::const_iterator2 it = iter.begin(); it != iter.end(); ++it) {
        Jprop(rows1,it.index2())=J(rows2,it.index2());
        // std::cout << "dp f_6 (" << rows1 << ", " << it.index2() << ") <- (" << rows2 << ", " << it.index2() << ") : 26 \n";
      }
    }
  }


  //* fdist7 */

  // right edge with upper right corner
  for(unsigned int y = n_y_ - 1; y >= 1;y--) {
    rows1 = GetPdfIndex(n_x_-1,y,Q_SW);
    rows2 = GetPdfIndex(n_x_-2,y-1,Q_SW);
    iter = J.find1(0, rows1, 0);

    for(mapped_matrix<double>::const_iterator2 it = iter.begin(); it != iter.end(); ++it) {
      Jprop(rows1,it.index2())= J(rows1,it.index2());
      // std::cout << "dp f_7 (" << rows1 << ", " << it.index2() << ") <- (" << rows1 << ", " << it.index2() << ") : 29 \n";
    }
    iter = J.find1(0, rows2, 0);

    for(mapped_matrix<double>::const_iterator2 it = iter.begin(); it != iter.end(); ++it) {
      Jprop(rows1,it.index2()) = Jprop(rows1,it.index2()) + J(rows2,it.index2());
      // std::cout << "dp f_7 (" << rows1 << ", " << it.index2() << ") <- (" << rows2 << ", " << it.index2() << ") : 30 \n";
    }
  }

  // upper edge without corners
  for(unsigned int x = n_x_ - 2; x >= 1; x--) {
    rows1 = GetPdfIndex(x,n_y_-1,Q_SW);
    rows2 = GetPdfIndex(x-1,n_y_-2,Q_SW);
    iter = J.find1(0, rows1, 0);

    for(mapped_matrix<double>::const_iterator2 it = iter.begin(); it != iter.end(); ++it) {
      Jprop(rows1,it.index2())= J(rows1,it.index2());
      // std::cout << "dp f_7 (" << rows1 << ", " << it.index2() << ") <- (" << rows1 << ", " << it.index2() << ") : 31 \n";
    }
    iter = J.find1(0, rows2, 0);

    for(mapped_matrix<double>::const_iterator2 it = iter.begin(); it != iter.end(); ++it) {
      Jprop(rows1,it.index2()) = Jprop(rows1,it.index2()) + J(rows2,it.index2());
      // std::cout << "dp f_7 (" << rows1 << ", " << it.index2() << ") <- (" << rows2 << ", " << it.index2() << ") : 32 \n";
    }
  }

  // inner elements
  for(unsigned int y = n_y_ - 2; y >= 1; y--) {
    for(unsigned int x= n_x_ - 2; x >= 1; x--) {
      rows1 = GetPdfIndex(x,y,Q_SW);
      rows2 = GetPdfIndex(x-1,y-1,Q_SW);
      iter = J.find1(0, rows2, 0);

      for(mapped_matrix<double>::const_iterator2 it = iter.begin(); it != iter.end(); ++it) {
        Jprop(rows1,it.index2())=J(rows2,it.index2());
        // std::cout << "dp f_7 (" << rows1 << ", " << it.index2() << ") <- (" << rows2 << ", " << it.index2() << ") : 33 \n";
      }
    }
  }


  //* fdist8 */

  // left edge with upper left corner
  for(unsigned int y = n_y_ - 1; y >= 1;y--) {
    rows1 = GetPdfIndex(0,y,Q_SE);
    rows2 = GetPdfIndex(1,y-1,Q_SE);
    iter = J.find1(0, rows1, 0);

    for(mapped_matrix<double>::const_iterator2 it = iter.begin(); it != iter.end(); ++it) {
      Jprop(rows1,it.index2())= J(rows1,it.index2());
      // std::cout << "dp f_8 (" << rows1 << ", " << it.index2() << ") <- (" << rows1 << ", " << it.index2() << ") : 36 \n";
    }
    iter = J.find1(0, rows2, 0);

    for(mapped_matrix<double>::const_iterator2 it = iter.begin(); it != iter.end(); ++it) {
      Jprop(rows1,it.index2()) = Jprop(rows1,it.index2()) + J(rows2,it.index2());
      // std::cout << "dp f_8 (" << rows1 << ", " << it.index2() << ") <- (" << rows2 << ", " << it.index2() << ") : 37 \n";
    }
  }

  // upper edge without corners
  for(unsigned int x = 1; x < n_x_-1; x++) {
    rows1 = GetPdfIndex(x,n_y_-1,Q_SE);
    rows2 = GetPdfIndex(x+1,n_y_-2,Q_SE);
    iter = J.find1(0, rows1, 0);

    for(mapped_matrix<double>::const_iterator2 it = iter.begin(); it != iter.end(); ++it) {
      Jprop(rows1,it.index2())= J(rows1,it.index2());
      // std::cout << "dp f_8 (" << rows1 << ", " << it.index2() << ") <- (" << rows1 << ", " << it.index2() << ") : 38 \n";
    }
    iter = J.find1(0, rows2, 0);

    for(mapped_matrix<double>::const_iterator2 it = iter.begin(); it != iter.end(); ++it) {
      Jprop(rows1,it.index2()) = Jprop(rows1,it.index2()) + J(rows2,it.index2());
      // std::cout << "dp f_8 (" << rows1 << ", " << it.index2() << ") <- (" << rows2 << ", " << it.index2() << ") : 39 \n";
    }
  }

  // inner elements
  for(unsigned int y = n_y_ - 2; y >= 1; y--) {
    for(unsigned int x = 1;x < n_x_-1; x++) {
      rows1 = GetPdfIndex(x,y,Q_SE);
      rows2 = GetPdfIndex(x+1,y-1,Q_SE);
      iter = J.find1(0, rows2, 0);

      for(mapped_matrix<double>::const_iterator2 it = iter.begin(); it != iter.end(); ++it) {
        Jprop(rows1,it.index2())=J(rows2,it.index2());
        // std::cout << "dp f_8 (" << rows1 << ", " << it.index2() << ") <- (" << rows2 << ", " << it.index2() << ") : 40 \n";
      }
    }
  }
}

Vector<double> LatticeBoltzmannPDE::d_pressuredrop_d_f(StdVector<double>& ux, StdVector<double>& uy, StdVector<double>& dloc)
{
  // Calculation of gradient of pressure drop with respect to design variable; By Georg Pingen, University of Colorado (CU), Boulder, Colorado
  double dA[9] = {0., 1., 0., -1., 0., 1., -1., -1., 1.};
  double dB[9] = {0., 0., 1., 0., -1., 1., 1., -1., -1.};
  double dUX[9],dUY[9];
  //boost::numeric::ublas::compressed_compressed_matrix<double> dPD(dpd->size1(),dpd->size2());
  // compressed_matrix<double> dPD(n_x_*n_y_*9,1,0.);
  mapped_matrix<double> dPD(n_elems * 9, 1, n_x_*n_y_*9);
  int index;
  for(unsigned int y=1;y<=n_y_;y++) {
    for(unsigned int x=1;x<=n_x_;x++) {
      index = (y-1)*n_x_+x-1;
      if(elements[index] == LBM_NODE_TYPE_INLET) // -2
      {
        for (int i=0;i<9;i++) {
          dUX[i]  = (dA[i] - ux[index])/dloc[index];
          dUY[i]  = (dB[i] - uy[index])/dloc[index];
          dPD(((y-1)*n_x_*9+(x-1)*9+i),0) = (1./3.+0.5*(ux[index]*ux[index] + uy[index]*uy[index])+dloc[index]*(ux[index]*dUX[i] + uy[index]*dUY[i])) / (double) inlet.GetSize();
        }
      }
      if(elements[index] == LBM_NODE_TYPE_OUTLET) // -3
      {
        for (int i=0;i<9;i++) {
          dUX[i]  = (dA[i] - ux[index])/dloc[index];
          dUY[i]  = (dB[i] - uy[index])/dloc[index];
          dPD(((y-1)*n_x_*9+(x-1)*9+i),0) = -(1./3.+0.5*(ux[index]*ux[index] + uy[index]*uy[index])+dloc[index]*(ux[index]*dUX[i] + uy[index]*dUY[i]))/static_cast<double>(outlet.GetSize());
        }
      }
    }
  }

  // compressed_matrix<double> dFdf(n_elems * 9,1,0.);
  mapped_matrix<double> dFdf(n_elems * 9, 1, n_elems * 9);

  d_propagate_d_f(dFdf,dPD);

  Vector<double> rhs(n_elems * 9);
  for(unsigned int i = 0, n = rhs.GetSize(); i < n; i++)
    rhs[i] = dFdf(i,0);

  return rhs; // no copy constructor
}

void LatticeBoltzmannPDE::SetPrecalculatedGradient(StdVector<DesignElement*>& design, Function* f)
{
  // the gradient over all lbm elements (including boundary) are computed as element wise scalar product of
  // x_parardiso * dRds where we could interpret dRds as vector but it is stored as sparse matrix and not all elements are contained


  // the solution of the adjoint system:
  Vector<double> x_pardiso2;
  x_pardiso2.Import("x_pardiso2.dat");

  LOG_DBG3(lbm_pde) << "SPG: x_pardiso2 -> " << x_pardiso2.ToString();

  // read dRds and store as vector, not set values keep 0.0
  Vector<double> dRds(n_q_ * n_elems, 0.0);
  // %%MatrixMarket matrix coordinate real symmetric
  // %
  // % Matrix exported by OLAS
  // %
  // 900 100 576
  // 100 12  0.000110019427034806
  // 101 12  -0.00388473301942947
  std::ifstream file("dRds.mtx");
  if(file.fail())
    throw Exception("cannot read dRds.mtx");

  std::string line;
  bool meta_info_read = false;
  while(std::getline(file, line))
  {
    if(line != "" && !boost::starts_with(line, "%"))
    {
      std::istringstream ss(line);
      // the first non-comment line is the meta data rows cols nnz
      if(!meta_info_read)
      {
        unsigned int rows, cols, nnz;
        if(!(ss >> rows >> cols >> nnz) || rows != n_q_*n_elems || cols != n_elems)
          throw Exception("error reading dRds.mtx, the first non-comment line has no proper meta data");
        meta_info_read = true;
      }
      else
      {
        unsigned int row, col;
        double val;
        if((ss >>  row >> col >> val))
          dRds[row - 1] = val;
        else
          throw Exception("error reading dRds.mtx values in out of line '" + line + "'");
      }
    }
  }
  file.close();

  LOG_DBG3(lbm_pde) << "SPG: dRds -> " << dRds.ToString();

  LOG_DBG3(lbm_pde) << "SPG: design size: " << design.GetSize();

  // set the actually requested gradients
  for(unsigned int e = 0; e < design.GetSize(); e++)
  {
    DesignElement* de = design[e];
    unsigned int idx = elem_to_idx[de->elem->elemNum]; // lbm idx
    double val = -1.0 * x_pardiso2.Inner(dRds, idx * n_q_, (idx + 1) * n_q_);
    de->AddGradient(f, val);
    LOG_DBG3(lbm_pde) << "SPG: idx=" << idx << " e=" << de->elem->elemNum << " xp=" << StdVector<double>::ToString(n_q_, x_pardiso2.GetPointer() + n_q_ * idx) << " dRds=" << StdVector<double>::ToString(n_q_, dRds.GetPointer() + n_q_ * idx) << " -> " << val;
  }
}


void LatticeBoltzmannPDE::InitTimeStepping()
{
  // timestepping formulation
  TS_alg_ = new PseudoTS( algsys_);
}

void LatticeBoltzmannPDE::CalcOutputCoupling() {

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

      if (quantity == LBM_VELOCITY) {
        solDeriv1_.SetAlgSysVector(getDeriv(FIRST_DERIV));
        solDeriv1_.NodeSolutionToCoupling((*values),*couplingnodes);
      }
      break;

    case ELEM:
      EXCEPTION("No Element coupling output");
    }

  }
}

bool LatticeBoltzmannPDE::HasOutput(SolutionType output) {

  if (output == LBM_VELOCITY || output == LBM_DENSITY)
    return true;
  return false;
}

void LatticeBoltzmannPDE::CalcResults( shared_ptr<BaseResult> res )
{
  // get the current state from the LBM Calculation.
  // we might come here AFTER LBM->Iterate() or during Iterate() when writing intermediate LBM iterations
  pdfs = lbm->GetPdfs();

  SolutionType solType = res->GetResultInfo()->resultType ;
  switch (solType) {
  case LBM_DENSITY:
    CalcDensities(res);
    break;
  case LBM_VELOCITY:
    CalcVelocities(res);
    break;
  case LBM_PRESSURE:
    CalcPressures(res);
    break;
  case MECH_PSEUDO_DENSITY:
  case LBM_PHYSICAL_PSEUDO_DENSITY:
    if(domain->GetErsatzMaterial(false) == NULL) // no exception
      res->Init();
    else
      domain->GetErsatzMaterial()->ExtractResults(res, false);
    break;
  case LBM_PROBABILITY_DISTRIBUTION:
    ExtractDistribution(res);
    break;

  default:
    WARN("Result type not computable by LatticeBoltzmann PDE" );
    break;
  }
}

double LatticeBoltzmannPDE::CalcLBMDensity(unsigned int idx) const
{
  double density = 0.0;
  for(unsigned int h = 0; h < n_q_; h++)
    density += pdf(idx,h);

  return density;
}

double LatticeBoltzmannPDE::CalcVelocityX(unsigned int idx, double density) const
{
  if (n_q_ == 9)
    return (pdf(idx, Q_E) + pdf(idx, Q_NE) + pdf(idx, Q_SE) - pdf(idx, Q_W) - pdf(idx, Q_NW) - pdf(idx, Q_SW)) / density;
  else
    return (pdf(idx, Q_NE) + pdf(idx, Q_E) + pdf(idx, Q_SE) + pdf(idx, Q_TE) + pdf(idx, Q_BE) - pdf(idx, Q_NW) - pdf(idx, Q_W) - pdf(idx, Q_SW) - pdf(idx, Q_TW) - pdf(idx, Q_BW)) / density;
}

double LatticeBoltzmannPDE::CalcVelocityY(unsigned int idx, double density) const
{
  if (n_q_ == 9)
    return (pdf(idx, Q_N)  + pdf(idx, Q_NE) + pdf(idx, Q_NW) - pdf(idx, Q_S) - pdf(idx, Q_SW) - pdf(idx, Q_SE)) / density;
  else
    return (pdf(idx, Q_N)  + pdf(idx, Q_NW) + pdf(idx, Q_NE) + pdf(idx, Q_TN) + pdf(idx, Q_BN) - pdf(idx, Q_S)- pdf(idx, Q_SE) - pdf(idx, Q_SW)  - pdf(idx, Q_TS) - pdf(idx, Q_BS)) / density;
}

double LatticeBoltzmannPDE::CalcVelocityZ(unsigned int idx, double density) const
{
  if (n_q_ == 9)
    return 0;
  else
    return (pdf(idx, Q_T) + pdf(idx, Q_TW) + pdf(idx, Q_TE) + pdf(idx, Q_TN) + pdf(idx, Q_TS) - pdf(idx, Q_B) - pdf(idx, Q_BW) - pdf(idx, Q_BE) - pdf(idx, Q_BN) - pdf(idx, Q_BS)) / density;
}

double LatticeBoltzmannPDE::CalcPressure(unsigned int idx) const
{
  double density = CalcLBMDensity(idx);
  double ux     = CalcVelocityX(idx, density);
  double uy     = CalcVelocityY(idx, density);
  double uz     = CalcVelocityZ(idx, density);

  return density / 3.0 + 0.5 * density * (ux * ux + uy * uy + uz * uz);
}


void LatticeBoltzmannPDE::CalcDensities( shared_ptr<BaseResult> base_result )
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


void LatticeBoltzmannPDE::CalcVelocities( shared_ptr<BaseResult> base_result )
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

    val[it.GetPos() * dim_]     = CalcVelocityX(idx, density);
    val[it.GetPos() * dim_ + 1] = CalcVelocityY(idx, density);
    if (dim_ == 3)
      val[it.GetPos() * dim_ + 2] = CalcVelocityZ(idx, density);
  }
}

void LatticeBoltzmannPDE::CalcPressures( shared_ptr<BaseResult> base_result )
{
Result<double>& res = dynamic_cast<Result<double>&>(*base_result);

  EntityIterator it = res.GetEntityList()->GetIterator();
  assert(it.GetType() == EntityList::ELEM_LIST);

  Vector<double>& val = res.GetVector();
  val.Resize(res.GetEntityList()->GetSize());

  // traverse the elements
  for(it.Begin(); !it.IsEnd(); it++)
    val[it.GetPos()] = CalcPressure(elem_to_idx[it.GetElem()->elemNum]);
}


double LatticeBoltzmannPDE::CalcPressureDrop()
{
  // Calculation of the pressure drop by Pingen
  double in = 0.0;
  for(unsigned int i = 0; i < inlet.GetSize(); i++)
    in += CalcPressure(inlet[i]);

  double out = 0.0;
  for(unsigned int i = 0; i < outlet.GetSize(); i++)
    out += CalcPressure(outlet[i]);

  return in / inlet.GetSize() - out / outlet.GetSize();
}

void LatticeBoltzmannPDE::ExtractDistribution( shared_ptr<BaseResult> base_result ){
  Result<double>& res = dynamic_cast<Result<double>&>(*base_result);

  EntityIterator it = res.GetEntityList()->GetIterator();
  assert(it.GetType() == EntityList::ELEM_LIST);
  Vector<double>& val = res.GetVector();
  val.Resize(res.GetEntityList()->GetSize() * n_q_);
  for(it.Begin(); !it.IsEnd(); it++)
  {
    for(unsigned int h = 0; h < n_q_; h++) {
      val[it.GetPos() * n_q_ + h] = pdf(elem_to_idx[it.GetElem()->elemNum],h);
    }
  }
}


// ***********************************************************************
//   Obtain information on desired output quantities from parameter file
// ***********************************************************************
void LatticeBoltzmannPDE::ReadSpecialResults()
{
}


void LatticeBoltzmannPDE::DefineAvailResults() {
  // =====================================================================
  // set solution information
  // =====================================================================
  PtrParamNode resultsNode = myParam_->Get("storeResults", ParamNode::PASS );

  // === solution by LBM solver ===
  shared_ptr<ResultInfo> prob(new ResultInfo);
  prob->resultType = LBM_PROBABILITY_DISTRIBUTION;
  StdVector<std::string> probDofNames;
  probDofNames.Resize(n_q_); // "f_0", "f_1", ....
  for(unsigned int i=0, n = n_q_; i < n; i++ )
    probDofNames[i] = "f_" + boost::lexical_cast<std::string>(i);
//  probDofNames[0] = "C";
//  probDofNames[1] = "E";
//  probDofNames[2] = "W";
//  probDofNames[3] = "N";
//  probDofNames[4] = "S";
//  probDofNames[5] = "T";
//  probDofNames[6] = "B";
//  probDofNames[7] = "NE";
//  probDofNames[8] = "SW";
//  probDofNames[9] = "NW";
//  probDofNames[10] = "SE";
//  probDofNames[11] = "TN";
//  probDofNames[12] = "BS";
//  probDofNames[13] = "TS";
//  probDofNames[14] = "BN";
//  probDofNames[15] = "TE";
//  probDofNames[16] = "BW";
//  probDofNames[17] = "TW";
//  probDofNames[18] = "BE";
  prob->dofNames = probDofNames;
  prob->unit = "";
  prob->entryType = ResultInfo::VECTOR;
  prob->fctType = shared_ptr<ConstFct>(new ConstFct());
  prob->definedOn = ResultInfo::ELEMENT;
  results_.Push_back(prob);
  availResults_.insert(prob);


  shared_ptr<ResultInfo> dens(new ResultInfo);
  dens->resultType = LBM_DENSITY;
  dens->unit =  "";
  dens->dofNames = "";
  dens->entryType = ResultInfo::SCALAR;
  dens->definedOn = ResultInfo::ELEMENT;
  dens->fctType = shared_ptr<ConstFct>(new ConstFct() );
  availResults_.insert( dens );


  shared_ptr<ResultInfo> pres(new ResultInfo);
  pres->resultType = LBM_PRESSURE;
  pres->unit =  "";
  pres->dofNames = "";
  pres->entryType = ResultInfo::SCALAR;
  pres->definedOn = ResultInfo::ELEMENT;
  pres->fctType = shared_ptr<ConstFct>(new ConstFct() );
  availResults_.insert(pres);


  shared_ptr<ResultInfo> velo(new ResultInfo);
  velo->resultType = LBM_VELOCITY;
  StdVector<std::string> velDofNames;
  if (dim_ == 2)
    velDofNames = "x", "y";
  else
    velDofNames = "x", "y", "z";
  velo->dofNames = velDofNames;
  velo->unit =  "";
  velo->entryType = ResultInfo::VECTOR;
  velo->definedOn = ResultInfo::ELEMENT;
  velo->fctType = shared_ptr<ConstFct>(new ConstFct() );
  availResults_.insert( velo );

  // === PSEUDO DENSITY for SIMP ===
  shared_ptr<ResultInfo> mechPD(new ResultInfo);
  mechPD->resultType = MECH_PSEUDO_DENSITY;
  mechPD->dofNames = "";
  mechPD->unit = "";
  mechPD->entryType = ResultInfo::SCALAR;
  mechPD->definedOn = ResultInfo::ELEMENT;
  mechPD->fctType = shared_ptr<ConstFct>(new ConstFct() );
  availResults_.insert( mechPD );

  shared_ptr<ResultInfo> physicalPD(new ResultInfo);
  physicalPD->resultType = LBM_PHYSICAL_PSEUDO_DENSITY;
  physicalPD->dofNames = "";
  physicalPD->unit = "";
  physicalPD->entryType = ResultInfo::SCALAR;
  physicalPD->definedOn = ResultInfo::ELEMENT;
  physicalPD->fctType = shared_ptr<ConstFct>(new ConstFct() );
  availResults_.insert( physicalPD );
}

// ***********************************************************************
//   Read LBM velocities from LBM output file
// ***********************************************************************
void LatticeBoltzmannPDE::ReadProbabilityDistribution(const std::string& filename)
{
  std::ifstream file(filename.c_str());
  if(file.fail())
    throw Exception("cannot open lbm result file " + std::string(filename));

  std::string line;
  unsigned int i = 0;
  while(std::getline(file, line) && i < n_elems)
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
    EXCEPTION("read " << i << " lines in " << filename << " but expected " << n_elems);


  file.close();

  // LOG_DBG3(lbm_pde) << "RD(" << filename << ") -> " << StdVector<double>::ToString(n_elems * 9, pdfs, 0);
}


void LatticeBoltzmannPDE::ExportCFS2LBM(const StdVector<double>& elements)
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

  for (UInt i = 0; i < inlet.GetSize(); ++i)
    file <<  u_x_ << " " << u_y_ << " " << u_z_ << std::endl;

  file.close();

  std::cout << "++ CFS2LBM.dat created" << std::endl;
}

void LatticeBoltzmannPDE::ExportMultipleFiles(const StdVector<double>& elements)
{
  std::ofstream por("por.dat");
  assert(n_elems == elements.GetSize());
  for(unsigned int i = 0, n = elements.GetSize(); i < n; i++)
    por << std::max(0.0, elements[i]) << std::endl; // boundary, inlet, outlet is all 0.0 in the old interface
  por.close();

  std::ofstream obst("obst.dat");
  // assume to have the mesh ordered from lower left starting to the right (lexicographic)
  // the OLD! interface assumes origin left upper and x downwards and y to the right!
  // TODO! Do by neighbours!
  Matrix<int> obst_m(n_y_, n_x_); // FIXME check order!
  for(unsigned int x = 0; x < n_x_; x++)
  {
    for(unsigned int y = 0; y < n_y_; y++)
    {
      // out: -1 bb, -2 inlet, -3 outlet, 0 ... 1 porisity
      // out: 1 bb1, 2 inlet, 3 outlet, 0 for porosity
      int v = (int) -1.0 * elements[n_x_ * y + x];
      obst << std::max(0, v) << " "; // org porosity would be -1
      obst_m(y, x) = std::max(0, v);  // origin of the matrix is left upper!
    }
    obst << std::endl;
  }
  obst.close();

  LOG_DBG2(lbm_pde) << "EMF: obst matrix (!= obst.dat):\n" << obst_m.ToString(0, true);

  std::ofstream ns("non_sing.dat");
  // add one to be one based
  for(unsigned int i = 0, n = non_sing.GetSize(); i < n; i++)
    ns << (non_sing[i] + 1) << " " << std::endl; // space to allow diff with matlab's non_sing.dat
  ns.close();

  std::ofstream data("data.dat");
  data << n_x_ << std::endl;
  data << n_y_ << std::endl;
  data << 1.0 << std::endl; // penalty parameter
  data << u_x_ << std::endl;
  data << u_y_ << std::endl;
  data << 1.0 << std::endl;  // density at the outlet
  data << omega_ << std::endl; // relaxation parameter for collision step
  data << maxIter_ << std::endl;
  data << convergence_ << std::endl; // lbm convergence tolerance
  data << non_sing.GetSize() << std::endl; // number of lines in non_sing.dat
  data << 1 << std::endl; // id of objective (1=pressure drop)
  data.close();


}

void LatticeBoltzmannPDE::WriteMatrix(const std::string& file, const compressed_matrix<double> & M)
{
  std::ofstream output(file.c_str());
  for (compressed_matrix<double>::const_iterator1 it1 = M.begin1(); it1 != M.end1(); ++it1) {
    for (compressed_matrix<double>::const_iterator2 it2 = it1.begin(); it2 != it1.end(); ++it2) {
      output << it2.index1()+1 << " " << it2.index2()+1 << "        " << *(it2) << std::endl;
    }
    output << std::endl;
  }
  output.close();
}

std::string LatticeBoltzmannPDE::ToString(const StdVector<double>& elements, bool x_fast, bool as_int) const
{
  std::stringstream ss;

  for(int y = n_y_ - 1; y >=  0; y--)
  {
    for(unsigned int x = 0; x < n_x_; x++)
    {
      int idx = x_fast ? y * n_x_  + x : x * n_y_ + y;
      if(as_int)
        ss << (int) std::abs(elements[idx]) << " ";
      else
        ss << std::abs(elements[idx]) << " ";
    }
    ss << " <- " << y << std::endl;
  }
  return ss.str();
}


void LatticeBoltzmannPDE::SetupElements()
{
  Grid* grd = domain->GetGrid();

  // find out which elements are inlet, outlet, boundary and inner ones
  StdVector<Elem*> elems;
  StdVector<Elem*> boundaries;
  // auxiliary vector
  // vector initialized with porosity value of inner cells
  elements.Resize(n_elems, 0.0);

  DesignSpace* space = domain->GetErsatzMaterial(false);
  if(space != NULL)
  {
    for(unsigned int r = 0; r < design_reg_.GetSize(); r++)
    {
      StdVector<Elem*> elems;
      grd->GetElems(elems, design_reg_[r]);

      for(unsigned int e = 0; e < elems.GetSize(); e++)
      {
        const Elem* elem = elems[e];
        int idx = space->Find(elem, true);
        double val = space->GetErsatzMaterialFactor(idx, Optimization::LBM);
        // the ordering of the design elements and the LBM elements might be different as it is defined for LBM
        // and might be arbitrary in the mesh which defines the ordering in the optimization design
        elements[elem_to_idx[elem->elemNum]] = val;
      }
    }
  }

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
  for (unsigned int i = 0; i < boundaries.GetSize(); ++i)
    elements[boundaries[i]->elemNum - 1] = -1.0;

  for(unsigned int i = 0; i < inlet.GetSize(); ++i)
    elements[inlet[i]] = -2.0;
  for(unsigned int i = 0; i < outlet.GetSize(); ++i)
    elements[outlet[i]] = -3.0;
}



void LatticeBoltzmannPDE::SetNonSingualrityIndices()
{
  assert(non_sing.IsEmpty());
  assert(!elements.IsEmpty());

  // reimplemented matlab code
  // we have the vector to export it as non_sing.dat for the old interface
  non_sing.Reserve(n_q_ * n_x_ * n_y_ * n_z_);



  for(unsigned int y = 0; y < n_y_; y++)
  {
    for(unsigned int x = 0; x < n_x_; x++)
    {
      unsigned int base = y * n_x_ * n_q_ + x * n_q_; // 2D
      //  LOG_DBG3(lbm_pde) << "INSI: test y=" << y << " x=" << x << " k=" << k << " base=" << base << " val -> " << val;
      if(elements[index(x,y)] != LBM_NODE_TYPE_BB) // no bounce back
      {
        // fluid node distribution functions are not deleted
        for(unsigned int h = 0; h < n_q_; h++) // 2D
          non_sing.Push_back(base + h);
      }
      else
      {
        // outpointing boundary distributions are not inserted and
        // distributions which point towards a bounce-back boundary point

        // test, if the node in direction f_1 (x+1) is not boundary but interior, inlet or outlet
        if(x+1 < n_x_ && elements[index(x+1, y)] != LBM_NODE_TYPE_BB) // if (inew<=lx && obst(inew,j)~=1)
          non_sing.Push_back(base + 1); // x(k)=(j-1)*lx*9+(i-1)*9+2;

        if(y+1 < n_y_ && elements[index(x,y+1)] != LBM_NODE_TYPE_BB) // y+1 = f_2_
          non_sing.Push_back(base + 2);

        if(x > 0 && elements[index(x-1,y)] != LBM_NODE_TYPE_BB) {
          non_sing.Push_back(base + 3);
        }

        if(y > 0 && elements[index(x,y-1)] != LBM_NODE_TYPE_BB)
          non_sing.Push_back(base + 4);

        if(x+1 < n_x_ && y+1 < n_y_ && elements[index(x+1,y+1)] != LBM_NODE_TYPE_BB)
          non_sing.Push_back(base + 5);

        if(x > 0 && y+1 < n_y_ && elements[index(x-1,y+1)] != LBM_NODE_TYPE_BB)
          non_sing.Push_back(base + 6);

        if(x > 0 && y > 0 && elements[index(x-1,y-1)] != LBM_NODE_TYPE_BB)
          non_sing.Push_back(base + 7);

        if(x+1 < n_x_ && y > 0 && elements[index(x+1,y-1)] != LBM_NODE_TYPE_BB)
          non_sing.Push_back(base + 8);
      }
    }
  }
}

void LatticeBoltzmannPDE::create_output(const char * file)
{
  // for debug purposes
    std::fstream f;
    f.precision(16);
    f.open(file, std::ios::out);
    for (unsigned int i=0;i< n_elems;i++)
    {
      for (unsigned int j=0;j<n_q_;j++) {
        f<<pdf(i,j)<<" ";
      }
      f<<std::endl;
    }
    f.close();
}

// void LatticeBoltzmannPDE::ToFile(const std::string& file, const compressed_matrix<double>& M)
void LatticeBoltzmannPDE::ToFile(const std::string& file, const mapped_matrix<double>& M)
{
  //writes a sparse matrix to file
  int number = M.nnz();
  double val;
  std::stringstream ss;
  ss.precision(16);
  // for(compressed_matrix<double>::const_iterator1 it = M.begin1(); it != M.end1(); ++it)
  for(mapped_matrix<double>::const_iterator1 it = M.begin1(); it != M.end1(); ++it)
  {
    // for(compressed_matrix<double>::const_iterator2 it2 = it.begin(); it2 != it.end(); ++it2)
    for(mapped_matrix<double>::const_iterator2 it2 = it.begin(); it2 != it.end(); ++it2)
    {
      val = M(it2.index1(),it2.index2());
      if(abs(val) > 1e-20)
        ss << (it2.index1() + 1)  << "\t" << (it2.index2() + 1) << "\t" << (val) << std::endl;
      else
        number--;
    }
  }
  std::fstream f;
  f.open(file.c_str(), std::ios::out);

  f << "%%MatrixMarket matrix coordinate real general\n";
  f << "%\n";
  f << "% Matrix extracted from boost::sparse_matrix\n";
  f << "%\n";
  f << M.size1() << "\t" << M.size2() << "\t" << number <<  "\n";
  f << ss.str() << std::endl;
  f.close();
}

