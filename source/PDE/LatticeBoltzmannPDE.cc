// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <stddef.h>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <cmath>
#include <filesystem>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "LatticeBoltzmannPDE.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ProgramOptions.hh"
#include "DataInOut/ResultHandler.hh"

#include "Domain/CoefFunction/CoefFunction.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "Domain/ElemMapping/Elem.hh"
#include "Domain/Domain.hh"
#include "Domain/ElemMapping/EntityLists.hh"
#include "Domain/Mesh/Grid.hh"
#include "Domain/Results/ResultInfo.hh"

#include "Utils/StdVector.hh"

#include "Driver/Assemble.hh"
#include "Driver/BaseDriver.hh"
#include "Driver/FormsContexts.hh"
#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Forms/BiLinForms/SingleEntryBiLinInt.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/Operators/GradientOperator.hh"

#include "FeBasis/BaseFE.hh"
#include "FeBasis/H1/FeSpaceH1Nodal.hh"
#include "General/Exception.hh"
#include "PDE/SinglePDE.hh"
#include "PDE/StdPDE.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/TransferFunction.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/ShapeMapDesign.hh"
#include "Optimization/Objective.hh"
#include "MatVec/Vector.hh"
#include "MatVec/StdMatrix.hh"
#include "MatVec/CRS_Matrix.hh"
#include "MatVec/SBM_Matrix.hh"
#include "MatVec/SBM_Vector.hh"

namespace CoupledField {

  using boost::numeric::ublas::matrix;


  class BaseMaterial;
  class SingleVector;

  // declare logging stream
  DEFINE_LOG(lbm_pde, "lbmpde")

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


  LatticeBoltzmannPDE::LatticeBoltzmannPDE(Grid* grid, PtrParamNode pn,
      PtrParamNode infoNode,
      shared_ptr<SimState> simState,
      Domain* domain ) : SinglePDE( grid, pn, infoNode,simState, domain )
  {
    forwardSim_.reset(new Timer("forwardSim"));
    backwardSim_.reset(new Timer("backwardSim"));
    setupAdjoint_.reset(new Timer("setupAdjoint"));
    setupPardiso_.reset(new Timer("setupPardiso"));
    solveLSE_.reset(new Timer("solveLSE"));
    dProp_.reset(new Timer("dprop",true));
    dColl_.reset(new Timer("dcoll",true));
    rhsSetup_.reset(new Timer("rhsSetup",true));
    delSing_ .reset(new Timer("delSing",true));
    dPressDrop_.reset(new Timer("dPressureDrop",true));

    pdename_ = "LatticeBoltzmann";
    pdematerialclass_ = MECHANIC;
    nonLin_ = false;
    lbm = NULL;
    numWriteResults_ = 0;
    numIterations_ = 0;

    method_ = "mechanic";

    // LBM parameters
    maxWallTime_ = myParam_->Get("LBM/maxWallTime")->As<double>();
    maxIter_     = myParam_->Get("LBM/maxIter")->As<unsigned int>();
    convergence_ = myParam_->Get("LBM/convergence")->As<double>();
    writeFrequency_ = myParam_->Get("LBM/writeFrequency")->As<unsigned int>();

    bool plot    = myParam_->Get("LBM/plot")->As<bool>();

    PtrParamNode bcsl = myParam_->Get("bcsAndLoads");
    // only one inlet region
    u_max_x_ = bcsl->HasByVal("inlet", "dof", "x") ? bcsl->GetByVal("inlet", "dof", "x")->Get("value")->As<double>() : 0.0;
    u_max_y_ = bcsl->HasByVal("inlet", "dof", "y") ? bcsl->GetByVal("inlet", "dof", "y")->Get("value")->As<double>() : 0.0;
    u_max_z_ = bcsl->HasByVal("inlet", "dof", "z") ? bcsl->GetByVal("inlet", "dof", "z")->Get("value")->As<double>() : 0.0;

    iface.SetName("LatticeBoltzmannPDE::Iface");
    iface.Add(INTERNAL, "internal");
    iface.Add(EXTERNAL, "external");
    iface_ = iface.Parse(pn->Get("LBM/solver")->As<std::string>());

    adjSRT_ = iface.Parse(pn->Get("LBM/solveAdjoint")->As<std::string>());

    InitRegions(pn, grid);

    // for 2D n_z_=1
    Vector<unsigned int> n = ShapeMapDesign::SetupLexicographicMesh(grid, boundary_reg_, elem_to_idx, idx_to_elem);
    n_x_ = n[0];
    n_y_ = n[1];
    n_z_ = n[2];
    n_elems = n_x_ * n_y_ * n_z_;

    // n_q gives number of discrete velocities of LBM model
    if (domain->GetGrid()->GetDim() == 2) {
      n_q_ = 9;
      dim_ = 2;
    }
    else if (domain->GetGrid()->GetDim() == 3) {
      n_q_ = 19;
      dim_ = 3;
    }

    StdVector<Elem*> tmp;
    grid->GetElemsByName(tmp,"inlet");
    inlet.Resize(tmp.GetSize());
    for(unsigned int i = 0; i < inlet.GetSize(); ++i)
      inlet[i] = elem_to_idx[tmp[i]->elemNum];

    // divide inlet velocity over number of inlet elements
    u_max_x_ /= inlet.GetSize();
    u_max_y_ /= inlet.GetSize();
    u_max_z_ /= inlet.GetSize();

    grid->GetElemsByName(tmp,"outlet");
    outlet.Resize(tmp.GetSize());
    for(unsigned int i = 0; i < outlet.GetSize(); ++i)
      outlet[i] = elem_to_idx[tmp[i]->elemNum];

    // if parabolic inflow profile required
    if (myParam_->Get("LBM/inflowProfile")->As<std::string>() == "parabolic") {
      parabolicInflow_ = true;

      SetupParabolicInflow();
    }
    else
      parabolicInflow_ = false;

    if (myParam_->Has("LBM/omega") && myParam_->Has("LBM/Re")) { // if Re and omega are given, adjust u_max_x
      omega_ = myParam_->Get("LBM/omega")->As<double>();
      Re_  = myParam_->Get("LBM/Re")->As<double>();

      u_max_x_ = Re_*(1/omega_ - 0.5) / (3.0 * inlet.GetSize());
      if (u_max_x_ / sqrt(3.0) > 0.2) {
        EXCEPTION("Mach number >= 0.2! Choose different Re or omega.");
      }
    }
    else if (myParam_->Has("LBM/omega")) { // if only omega given, use given inlet vel. and omega to calculate Re
      omega_       = myParam_->Get("LBM/omega")->As<double>();
      //calculate Reynolds number of fluid flow
      // re = inletSize * |u| / kin. viscosity
      double u_mean_x = u_max_x_;
      double u_mean_y = u_max_y_;
      double u_mean_z = u_max_z_;
      if (parabolicInflow_) {
        u_mean_x = 2.0/3.0 * u_max_x_;
        u_mean_y = 2.0/3.0 * u_max_y_;
        u_mean_z = 2.0/3.0 * u_max_z_;
      }

      Re_ = inlet.GetSize() * sqrt(u_mean_x * u_mean_x + u_mean_y * u_mean_y+ u_mean_z * u_mean_z) / (1/3.0 * (1/omega_ - 0.5));
    }
    else if (myParam_->Has("LBM/Re")) { // if only Re given, calc omega
      if (parabolicInflow_)
        omega_ = 1.0 / ( 3*inlet.GetSize() * sqrt(1.5*1.5*(u_max_x_ * u_max_x_ + u_max_y_ * u_max_y_+ u_max_z_ * u_max_z_)) / Re_ + 0.5);
      else
        omega_ = 1.0 / ( 3*inlet.GetSize() * sqrt(u_max_x_ * u_max_x_ + u_max_y_ * u_max_y_+ u_max_z_ * u_max_z_) / Re_ + 0.5);
    }

    if (myParam_->Has("LBM/omega")) {
      omega_       = myParam_->Get("LBM/omega")->As<double>();
      //calculate Reynolds number of fluid flow
      // re = inletSize * |u| / kin. viscosity
      double u_mean_x = u_max_x_;
      double u_mean_y = u_max_y_;
      double u_mean_z = u_max_z_;
      if (parabolicInflow_) {
        u_mean_x = 2.0/3.0 * u_max_x_;
        u_mean_y = 2.0/3.0 * u_max_y_;
        u_mean_z = 2.0/3.0 * u_max_z_;
      }

      Re_ = inlet.GetSize() * sqrt(u_mean_x * u_mean_x + u_mean_y * u_mean_y+ u_mean_z * u_mean_z) / (1/3.0 * (1/omega_ - 0.5));
    }
    else if (myParam_->Has("LBM/Re")) {
      Re_       = myParam_->Get("LBM/Re")->As<double>();
      if (parabolicInflow_)
        omega_ = 1.0 / ( 3*inlet.GetSize() * sqrt(1.5*1.5*(u_max_x_ * u_max_x_ + u_max_y_ * u_max_y_+ u_max_z_ * u_max_z_)) / Re_ + 0.5);
      else
        omega_ = 1.0 / ( 3*inlet.GetSize() * sqrt(u_max_x_ * u_max_x_ + u_max_y_ * u_max_y_+ u_max_z_ * u_max_z_) / Re_ + 0.5);
      if (omega_ >= 2)
        EXCEPTION("Omega=" << omega_ << " must be smaller 2. Choose different Reynolds number or inlet velocity!");
    }

    //Initializing storage for PDFs
    pdfs_.Resize(n_elems * n_q_);
    adjPdfs.Resize(n_elems * n_q_);
    adjPdfs.Init(0.0);

    // vector storing adjoint SRT collision matrices
    adjSRTCollision.Resize(n_elems);
    d_pdrop_d_f.Resize(n_elems);
    adjCollisions.Resize(n_elems);
    StdVector<double> zeros(n_q_);
    zeros.Init();
    // Initialize d_pdrop_d_f with zero vectors
    for (unsigned int elem = 0; elem < n_elems; elem++) {
      d_pdrop_d_f[elem] = zeros;
      adjCollisions[elem].Resize(n_q_*n_q_);
    }

    if(iface_ == INTERNAL) {
      lbm = new LatticeBoltzmann(dim_, n_x_, n_y_, n_z_, u_max_x_, u_max_y_, u_max_z_, u_in_, omega_, maxIter_, convergence_, plot, writeFrequency_);
    }

    microVelDirections_ = lbm->GetPDFDirectionVectors();
    invPDFDirections_ = lbm->GetinvPDFDirections();

//    PtrParamNode sensAnalysis = infoNode_->Get(ParamNode::SUMMARY)->Get("sensitvityAnalysis", ParamNode::INSERT);
//    sensAnalysis->Get("timer", ParamNode::APPEND)->SetValue(backwardSim_);
//    sensAnalysis->Get("timer", ParamNode::APPEND)->SetValue(setupAdjoint_);
//    sensAnalysis->Get("timer", ParamNode::APPEND)->SetValue(setupPardiso_);
//    sensAnalysis->Get("timer", ParamNode::APPEND)->SetValue(solveLSE_);
//    sensAnalysis->Get("timer", ParamNode::APPEND)->SetValue(propSetup_);
//    sensAnalysis->Get("timer", ParamNode::APPEND)->SetValue(rhsSetup_);
//    sensAnalysis->Get("timer", ParamNode::APPEND)->SetValue(delSing_);
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
    obstacle_reg_ = -1;
    design_reg_.Reserve(regions.GetSize() - 2); // ommit boundary region and inclusion region

    for(unsigned int i = 0; i < regions.GetSize(); i++)
    {
      RegionIdType reg = grid->GetRegion().Parse(regions[i]->Get("name")->As<std::string>());
      if(regions[i]->Get("boundary")->As<bool>())
      {
        if(boundary_reg_ != -1)
          throw Exception("only a single region might have the boundary attribute set");
        else
          boundary_reg_ = reg;
      }
      else if(regions[i]->Get("name")->As<std::string>() == "obstacle")
      {
        obstacle_reg_ = reg;
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
    // To be good this needs to handled by element neighbors!
    if(grid->GetNumElems() != n_elems)
      EXCEPTION("the current implementation assumes the whole mesh to be used for LBM and lexicographic ordering. Mesh has " << grid->GetNumElems() << " but we assume " << n_elems);

    idx_to_elem.Resize(n_elems);
    for(unsigned int i = 0; i < n_elems; i++)
      idx_to_elem[i] = i+1; // one based

    elem_to_idx.Resize(n_elems + 1); // one-based elem_nr
    for(unsigned int i = 0, n = elem_to_idx.GetSize(); i < n; i++)
      elem_to_idx[i] = i-1; // -1 for not appropriate idx
  }

  void LatticeBoltzmannPDE::DefineIntegrators()
  {
    ////////////////////// code taken from testPDE

    // help variables for parameter checking
    RegionIdType actRegion;
    BaseMaterial * actSDMat = NULL;

    // Define integrators for "standard" materials
    std::map<RegionIdType, BaseMaterial*>::iterator it;

    //get FEFunction and space
    shared_ptr<BaseFeFunction> feFunc = feFunctions_[TEST_DOF];
    shared_ptr<FeSpace> mySpace = feFunc->GetFeSpace();

//    for ( it = materials_.begin(); it != materials_.end(); it++ ) {
      // Set current region and material
      actRegion = materials_.begin()->first;
      actSDMat = materials_.begin()->second;

      // create new entity list
      // we use colList in order to be able to use constant FESpace, which is described by only one equation
      shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
      actSDList->SetRegion( actRegion );

      // get current region name and get grip of paramNode
      std::string actRegionName = ptGrid_->GetRegion().ToString( actRegion );

      // --- Set the FE ansatz for the current region ---
      PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",actRegionName.c_str());
      mySpace->SetRegionApproximation(actRegion, "default","default");

      // pass entitylist of fespace / fefunction
      feFunc->AddEntityList( actSDList );

      // ====================================================================
      // stiffness integrator (taken from TestPDE)
      // ====================================================================
      PtrCoefFct beta = actSDMat->GetScalCoefFnc( DENSITY, Global::REAL );

      SingleEntryBiLinInt* stiffInt = NULL;

      shared_ptr<CoefFunction> coefFunc(new CoefFunctionLBM<Double>(this,feFunc,results_[0]));
      stiffInt = new SingleEntryBiLinInt(n_q_,coefFunc, false); // we want non-symmetric system matrix
      stiffInt->SetName("StiffnessIntegrator");
      LOG_TRACE(lbm_pde) << "Integrator symmetric? " << stiffInt->IsSymmetric();

      BiLinFormContext * stiffIntDescr = new BiLinFormContext(stiffInt, STIFFNESS );

      stiffIntDescr->SetEntities( actSDList, actSDList );
      stiffIntDescr->SetFeFunctions(feFunc,feFunc);
      stiffInt->SetFeSpace( mySpace );

      assemble_->AddBiLinearForm( stiffIntDescr );
//    }
      PtrParamNode sensAnalysis = infoNode_->Get(ParamNode::SUMMARY)->Get("sensitvityAnalysis", ParamNode::INSERT);
      sensAnalysis->Get("timer", ParamNode::APPEND)->SetValue(forwardSim_);
      sensAnalysis->Get("timer", ParamNode::APPEND)->SetValue(backwardSim_);
      sensAnalysis->Get("timer", ParamNode::APPEND)->SetValue(setupAdjoint_);
      sensAnalysis->Get("timer", ParamNode::APPEND)->SetValue(setupPardiso_);
      sensAnalysis->Get("timer", ParamNode::APPEND)->SetValue(solveLSE_);
      sensAnalysis->Get("timer", ParamNode::APPEND)->SetValue(dProp_);
      sensAnalysis->Get("timer", ParamNode::APPEND)->SetValue(dColl_);
      sensAnalysis->Get("timer", ParamNode::APPEND)->SetValue(dPressDrop_);
      sensAnalysis->Get("timer", ParamNode::APPEND)->SetValue(rhsSetup_);
      sensAnalysis->Get("timer", ParamNode::APPEND)->SetValue(delSing_);
  }

  void LatticeBoltzmannPDE::InitCoupling(PDECoupling * coupling)
  {
  }


  std::map<SolutionType, shared_ptr<FeSpace> > LatticeBoltzmannPDE::CreateFeSpaces( const std::string&  formulation, PtrParamNode infoNode )
  {
    std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;
    if( formulation == "default"){
      PtrParamNode potSpaceNode = infoNode->Get("testDof"); // dummy state solution TEST_DOF
      crSpaces[TEST_DOF] = FeSpace::CreateInstance(myParam_,potSpaceNode,FeSpace::CONSTANT, ptGrid_);
      crSpaces[TEST_DOF]->Init(solStrat_);
    }else{
      EXCEPTION("The formulation " << formulation << "of the lattice Boltzmann PDE is not known!");
    }
    return crSpaces;
  }

  void LatticeBoltzmannPDE::DefineSolveStep()
  {
    solveStep_ = new StdSolveStep(*this);
  }

  void LatticeBoltzmannPDE::Solve()
  {
    // infoNode_ is not set yet in the constructor
    PtrParamNode in = infoNode_->Get(ParamNode::HEADER)->Get("LBM");
    in->Get("omega_nu")->SetValue(omega_);
    in->Get("Re")->SetValue(Re_);
    in->Get("maxIter")->SetValue(maxIter_);
    in->Get("maxWallTime")->SetValue(maxWallTime_);
    in->Get("convergence")->SetValue(convergence_);
    in->Get("iface")->SetValue(iface.ToString(iface_));
    in->Get("u_max_x")->SetValue(u_max_x_);
    in->Get("u_max_y")->SetValue(u_max_y_);
    in->Get("u_max_z")->SetValue(u_max_z_);
    // we divide the input velocity by inlet.GetSize()
    in->Get("u_max_scaling")->SetValue(inlet.GetSize());

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
      case EXTERNAL:
      {

        executable = myParam_->Get("LBM")->Get("lbm")->As<std::string>();

        ExportCFS2LBM(elements);

        if(!std::filesystem::exists(executable))
          EXCEPTION("Could not find executable '" + executable + "', might be not in path");

        std::cout << "++ Calling external LBM solver .. \n" << std::endl;
        std::string command = executable + " -o LBM2CFS.dat CFS2LBM.dat";
        int err = system(command.c_str());
        if (err)
          EXCEPTION("LBM simulation failed, no outputs available! \n");

        ReadProbabilityDistribution("LBM2CFS.dat");

        break;
      }
      case INTERNAL:
      {
        LOG_DBG(lbm_pde) << "call internal LBM";
        LOG_DBG2(lbm_pde) << "elements\n" << ToString(elements, true, true);

        // this data is the final simulation but it lacks an additional propagation step for the adjoint calculation
        // We need to perform an additional propagation step as base for the adjoint system setup
        forwardSim_->Start();
        StdVector<double>* tmp = lbm->Iterate(elements, in->Get("LBM"));

        pdfs_ = *tmp;
        in->Get("original_pressure_drop")->SetValue(CalcPressureDrop());

        pdfs_ = *tmp;
        lbm->Prop_step();

        pdfs_ = *tmp;
        in->Get("prop_step_pressure_drop")->SetValue(CalcPressureDrop());
        forwardSim_->Stop();

        numWriteResults_ = lbm->GetNumWriteResults();

        numIterations_ = lbm->GetNumIterations();
        break;
      }
    }

    timer.Stop();
    state_.Stop();

    in = infoNode_->Get(ParamNode::SUMMARY)->Get("stateProblem");
    in->Get("original_pressure_drop")->SetValue(CalcPressureDrop());
    in->Get("prop_step_pressure_drop")->SetValue(CalcPressureDrop());
    in->Get("totalTimer/cpu")->SetValue(state_.GetCPUTime());
    in->Get("totalTimer/wall")->SetValue(state_.GetWallTime());
    in->Get("totalTimer/calls")->SetValue(state_.GetCalls());

  }

  void LatticeBoltzmannPDE::SetupSensitivityAnalysis(StdVector<double>& ux, StdVector<double>& uy, StdVector<double>& uz, StdVector<double>& dloc, StdVector<double>& weights)
  {
    ux.Resize(n_elems);
    uy.Resize(n_elems);
    uz.Resize(n_elems);
    dloc.Resize(n_elems);

    // macroscopic values for each m_node
    for(unsigned int i = 0; i < n_elems; i++)
    {
      double density = CalcLBMDensity(i);
      ux[i] = CalcVelocityX(i, density);
      uy[i] = CalcVelocityY(i, density);
      uz[i] = CalcVelocityZ(i, density);
      dloc[i] = density;

      if (dim_ == 2)
        assert(uz[i] == 0);

      // sets quadrature weights
      weights.Resize(n_q_);
      if (dim_ == 3) {
        weights[Q_0] = 1. / 3.;
        weights[Q_E] = 1. / 18.;
        weights[Q_N] = 1. / 18.;
        weights[Q_W] = 1. / 18.;
        weights[Q_S] = 1. / 18.;
        weights[Q_T] = 1. / 18.;
        weights[Q_B] = 1. / 18.;
        weights[Q_NE] = 1. / 36.;
        weights[Q_NW] = 1. / 36.;
        weights[Q_SW] = 1. / 36.;
        weights[Q_SE] = 1. / 36.;
        weights[Q_TN] = 1. / 36.;
        weights[Q_BS] = 1. / 36.;
        weights[Q_TS] = 1. / 36.;
        weights[Q_BN] = 1. / 36.;
        weights[Q_TE] = 1. / 36.;
        weights[Q_BW] = 1. / 36.;
        weights[Q_TW] = 1. / 36.;
        weights[Q_BE] = 1. / 36.;
      }
      else {
        weights[Q_0] = 4.0 / 9.0;
        weights[Q_E] = 1.0 /  9.0;
        weights[Q_N] = 1.0 /  9.0;
        weights[Q_W] = 1.0 /  9.0;
        weights[Q_S] = 1.0 /  9.0;
        weights[Q_NE] = 1.0 / 36.0;
        weights[Q_NW] = 1.0 / 36.0;
        weights[Q_SW] = 1.0 / 36.0;
        weights[Q_SE] = 1.0 / 36.0;
      }
    }
  }


void LatticeBoltzmannPDE::SensitivityAnalysis(TransferFunction* tf, Function* f, DesignSpace* space)
{
  // needs to be called after Solve() !!
  setupAdjoint_->Start();
  dColl_->Start();

  // Initialization of the data structures
  StdVector<double> ux(n_elems);
  StdVector<double> uy(n_elems);
  StdVector<double> uz(n_elems);
  StdVector<double> dloc(n_elems);
  StdVector<double> weights(n_q_);

  // setting values in data structures
  SetupSensitivityAnalysis(ux, uy, uz, dloc, weights);

  // derivative of the residual with respect to design variables
  Vector<double> dRds(n_q_ * n_elems);
  dRds.Init(0.0);

  // derivative of the collision operator with respect to p
  mapped_matrix<double> col_jacobi(n_q_ * n_elems, n_q_ * n_elems);

  StdVector<double> dfeqdux(n_q_);
  StdVector<double> dfeqduy(n_q_);
  StdVector<double> dfeqduz(n_q_);

  Matrix<double> block(n_q_,n_q_);

  for(unsigned int index = 0; index < n_elems ; index++)
  {
    block.InitValue(0.0);
    // interior
    if(elements[index] >= 0 && !obstacles.Contains(index))  {
      // Jacobi matrix of the collision operator with respect to the design variables
      d_collision_step_d_f(index, block, dfeqdux, dfeqduy, dfeqduz, ux, uy, uz, dloc, weights);

      // simple choice of the mapping between p and s
      DesignElement* de = space->Find(idx_to_elem[index], DesignElement::DENSITY);
      double scale = -1.0 * tf->Derivative(de, DesignElement::SMART);
      for (unsigned int k = 0; k < n_q_; k++)
        dRds[GetPdfIndex(index,k)] = omega_ * (dfeqdux[k] * scale * ux[index] + dfeqduy[k] * scale * uy[index] + dfeqduz[k] * scale * uz[index]);
    }
    else if(elements[index] == LBM_NODE_TYPE_BB) {
      d_bounceback_d_f(index,block); // bounce-back sensitivities
    }
    else if(elements[index] == LBM_NODE_TYPE_INLET) {
      d_inflow_d_f(index,block, weights); // derivative at inlet with respect to f
    }
    else if(elements[index] == LBM_NODE_TYPE_OUTLET) {
      d_outflow_d_f(index, block, ux, uy, uz, dloc, weights); // derivative at outlet with respect to f
    }
    else if (elements[index] != LBM_NODE_TYPE_OBSTACLE){
      //      assert(false);
    }

    // fill transpose of block in col_jacobi
    for(unsigned int k = 0; k < n_q_; k++)
      for(unsigned int l = 0; l< n_q_; l++) {
        col_jacobi(GetPdfIndex(index,k),GetPdfIndex(index,l)) = block[l][k];
      }
  }

  dColl_->Stop();
  dPressDrop_->Start();

  // right-hand side
  Vector<Double> b = d_pressuredrop_d_f(ux, uy, uz);

  dPressDrop_->Stop();

  // using external LSE solver
  if (adjSRT_ == EXTERNAL)
  {
    dProp_->Start();
    // the real jacobi combines d_collision with d_propagation
    mapped_matrix<double> Jacobi(n_q_ * n_elems , n_q_ * n_elems);
    d_propagate_d_f(Jacobi,col_jacobi);

    // substract identity matrix from LBM derivatives: dR/dr = d(LBM)/df - 1
    for(unsigned int i=0;i < n_q_ * n_elems; i++)
      Jacobi(i,i) -= 1.;

    dProp_->Stop();

    delSing_->Start();
    //Data for Pardiso solver
    compressed_matrix<double> Jacobi_new(n_elems * n_q_, n_elems * n_q_, Jacobi.nnz());

    // Delete singular rows from the Jacobian
    DeleteSingularities(Jacobi,Jacobi_new);

    delSing_->Stop();
    setupAdjoint_->Stop();

    setupPardiso_->Start();
    // use the cfs system matrix to solve the adjoint system
    StdMatrix* mat = algsys_->GetMatrix(SYSTEM)->GetPointer(0,0);
    LOG_DBG(lbm_pde) << "SA: " << mat->ToString(',','\n');
    LOG_DBG(lbm_pde) << "SA: size of system's matrix=" << mat->GetNumRows() << " x " << mat->GetNumCols();
    LOG_DBG(lbm_pde) << "SA: mat structure=" << mat->GetStructureType();
    LOG_DBG(lbm_pde) << "SA: mat storage=" << mat->GetStorageType();
    LOG_DBG(lbm_pde) << "SA: mat entry type=" << mat->GetEntryType();

    CRS_Matrix<double>* crs = dynamic_cast<CRS_Matrix<double>*>(mat);
    assert(crs != NULL);

    crs->SetSize(n_elems * n_q_, n_elems * n_q_, Jacobi_new.nnz());
    matrix_sparse_to_crs(Jacobi_new, crs->GetDataPointer(), crs->GetRowPointer(), crs->GetColPointer());

//    // time to setup adjoint system before solving
//    double setup_wall = timer.GetWallTime();
//    double setup_cpu = timer.GetCPUTime();

    LOG_DBG(lbm_pde) << "SA: d_pressuredrop_d_f=" << b.ToString();
    SBM_Vector rhs(BaseMatrix::DOUBLE);
    rhs.Resize(1);
    rhs.SetSubVector(GenerateSingleVectorObject(BaseMatrix::DOUBLE,b.GetSize()),0);
    rhs.AddToSubVector(b,0);

    LOG_TRACE(lbm_pde) << "SA: " << " size of rhs: " << rhs.GetPointer(0)->GetSize();
    LOG_DBG3(lbm_pde) << "SA: " << " rhs=" << rhs.GetPointer(0)->ToString();

    algsys_->InitRHS(rhs);
    algsys_->SetupSolver();
    setupPardiso_->Stop();

    solveLSE_->Start();
    algsys_->Solve();
    Vector<double> sol;
    algsys_->GetSolutionVal(sol,0,false);
    solveLSE_->Stop();

    for(unsigned int e = 0; e < f->elements.GetSize(); e++)
    {
      DesignElement* de = f->elements[e];
      unsigned int idx = elem_to_idx[de->elem->elemNum]; // lbm idx
      double val = -1.0 * sol.Inner(dRds, idx * n_q_, (idx + 1) * n_q_);
      de->AddGradient(f, val);
    }

//    adjoint->Get("timer/cpu")->SetValue(timer.GetCPUTime());
//    adjoint->Get("timer/wall")->SetValue(timer.GetWallTime());
//    adjoint->Get("setupTimer/cpu")->SetValue(setup_cpu);
//    adjoint->Get("setupTimer/wall")->SetValue(setup_wall);
//    adjoint->Get("setupTimer/d_coll")->SetValue(d_collision_setup);
//    adjoint->Get("setupTimer/d_prop")->SetValue(d_propagation_setup - d_collision_setup);
//    adjoint->Get("setupTimer/del_sing")->SetValue(rhs_setup - delete_sing_setup);
//    adjoint->Get("setupTimer/rhs")->SetValue(delete_sing_setup - d_propagation_setup);
  }
  else // solving adjoint system with same complexity as LBM simulation
  {
    setupAdjoint_->Stop();
    backwardSim_->Start();
//    StdVector<double>* tmp = lbm->IterateAdjointSRT(infoNode_->Get(ParamNode::PROCESS)->Get("stateProblem/LBM/adjoint"),adjSRTCollision,d_pdrop_d_f);
    StdVector<double>* tmp = lbm->IterateAdjointSRT(infoNode_->Get(ParamNode::PROCESS)->Get("stateProblem/LBM/adjoint"),adjCollisions,d_pdrop_d_f);
    adjPdfs = *tmp;
    backwardSim_->Stop();

    for(unsigned int e = 0; e < f->elements.GetSize(); e++)
    {
      DesignElement* de = f->elements[e];
      unsigned int idx = elem_to_idx[de->elem->elemNum]; // lbm idx
      Vector<double> sol(n_q_), d_coll_d_s(n_q_);
      for (unsigned int dir = 0; dir < n_q_; dir++) {
        sol[dir] = adjPdfs[GetPdfIndex(idx,dir)];
        d_coll_d_s[dir] = dRds[GetPdfIndex(idx,dir)];
      }
      double val = -d_coll_d_s.Inner(sol);
      de->AddGradient(f, val);
    }
//    adjoint->Get("iterative/timer/cpu")->SetValue(timer.GetCPUTime());
//    adjoint->Get("iterative/timer/wall")->SetValue(timer.GetWallTime());
//    adjoint->Get("setupTimer/d_coll")->SetValue(d_collision_setup);

  }
//  timer.Stop();
//  adjoint_.Stop();
//  adjoint = infoNode_->Get(ParamNode::SUMMARY)->Get("adjoint");
//  adjoint->Get("totalTimer/cpu")->SetValue(adjoint_.GetCPUTime());
//  adjoint->Get("totalTimer/wall")->SetValue(adjoint_.GetWallTime());
//  adjoint->Get("totalTimer/calls")->SetValue(adjoint_.GetCalls());
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

/**
 * This function "deletes" the entries from the Jacobian causing its singularity. This is done by copying the valid values from the old Jacobian M to a new one.
 * - non_sing is a vector containing indices of non singular elements
 * - non_sing is being copied to temporary set nss due to performance reasons
 *
 * Basic procedure:
 * - Iterate over all rows of M
 * - If already row index is not contained in non_sing, we replace this row with the unit row vector to assure full rank of new matrix
 * - If row and column indices of M is contained in non_sing, this values is copied to new Jacobian
 */

void LatticeBoltzmannPDE::DeleteSingularities(const mapped_matrix<double> & M, compressed_matrix<double> & output)
{
  std::set<unsigned int> nss(non_sing.Begin(), non_sing.End()); // non_sing set, searching a set is only O(log(n)) instead of O(n) for a vector

  // delete singular rows of the Jacobian in the sensitivity analysis for the optimization
  double val;

  for(mapped_matrix<double>::const_iterator1 it = M.begin1(); it != M.end1(); ++it) { // iterate over all rows of matrix M

    if (nss.find(it.index1()) == nss.end()) { // if this row causes singularity of Jacobian
      // to assure full rank of Jacobian, fill rows that cause singularity with unit vector: diagonal entry gets value 1
      output(it.index1(),it.index1()) = 1.0;
      continue;
    }
    else { // if row index is in non_sing vector
      for(mapped_matrix<double>::const_iterator2 it2 = it.begin(); it2 != it.end(); ++it2) {
        val = M(it2.index1(),it2.index2());
        if (nss.find(it2.index2()) != nss.end()) { // if column index is also in non_sing vector
          //it2.index1()  row index , it2.index2() column index
          output(it2.index1(),it2.index2()) = val;
          LOG_DBG3(lbm_pde) << "DS: o("<< it2.index1() << ", " << it2.index2() << ") = " << val;
        }
      }
    }
  }
}


void LatticeBoltzmannPDE::d_collision_step_d_f(unsigned int index, Matrix<double>& block, StdVector<double>& dfeqdux, StdVector<double>& dfeqduy, StdVector<double>& dfeqduz, const StdVector<double>& ux, const StdVector<double>& uy, const StdVector<double>& uz, const StdVector<double>& dloc, const StdVector<double>& weight)
{
  // gradient of the collision step with respect to the design variables
  // partial derivative of f^eq with respect to rho, ux, uy and uz: CORRECT (FDM)
  StdVector<double> dfeqdrho(n_q_);
  double scale = 1. - elements[index]; // 1.-pow(por[0][index],penal);
  assert(scale >= 0);
  LOG_DBG3(lbm_pde) << "d_collision_step_d_f: index = " << index << " scale=" << scale;
  //LOG_DBG3(lbm_pde) << "d_collision_step_d_f: index = " << index << " pdf=" << StdVector<double>::ToString(9, &pdfs.GetPointer()[index * 9]);
  //partial derivative of f^eq with respect to rho, ux, uy and uz: CORRECT (FDM)
  double us1, us2, us3, dot, norm;
  LatticeBoltzmann::PDFDirectionVector transform;

  //gradient of u_x with respect to f
  StdVector<double> duxdf(n_q_);
  //gradient of u_y with respect to f
  StdVector<double> duydf(n_q_);
  //gradient of u_z with respect to f
  StdVector<double> duzdf(n_q_);

  double invdloc = 1.0 / dloc[index];

  for (unsigned int i = 0; i < n_q_; i++)
  {
    transform = microVelDirections_[i];
    us1 = scale * ux[index];
    us2 = scale * uy[index];
    us3 = scale * uz[index];

    dot = transform.off_x * us1 + transform.off_y * us2 + transform.off_z * us3;
    norm = us1 * us1 + us2 * us2 + us3 * us3;

    dfeqdrho[i] = weight[i] * (1. + 3.*dot + 4.5*dot*dot - 1.5* norm);
    dfeqdux[i] = weight[i] * dloc[index]*(3. * transform.off_x + 9. * transform.off_x * dot - 3. * us1);
    dfeqduy[i] = weight[i] * dloc[index]*(3. * transform.off_y + 9. * transform.off_y * dot - 3. * us2);
    dfeqduz[i] = weight[i] * dloc[index]*(3. * transform.off_z + 9. * transform.off_z * dot - 3. * us3);

    LOG_DBG3(lbm_pde) << "d_collision_step_d_f: w = " << weight[i] << " dloc=" << dloc[index] << " dot=" << dot;

    duxdf[i] = scale * invdloc * (-ux[index] + transform.off_x);
    duydf[i] = scale * invdloc * (-uy[index] + transform.off_y);
    duzdf[i] = scale * invdloc * (-uz[index] + transform.off_z);
    // automatic testing
    if (dim_ == 2) {
      assert(us3 == 0);
      assert(dot == transform.off_x * us1 + transform.off_y * us2);
      assert(norm == us1 * us1 + us2 * us2);
      assert(dfeqduz[i] == 0);
      assert(duzdf[i] == 0);
    }
  }

  // partial derivatives of f^eq with respect to f: CORRECT (FDM)
  StdVector< StdVector<double> > dfeqdf(n_q_);
  for (unsigned int dir = 0; dir < n_q_; dir++ ){
    dfeqdf[dir].Resize(n_q_);
  }

  for (unsigned int j = 0; j < n_q_; j++)
    for (unsigned int i = 0; i < n_q_; i++)
      dfeqdf[i][j] = dfeqdrho[i] + dfeqdux[i]*duxdf[j] + dfeqduy[i]*duydf[j] + dfeqduz[i] * duzdf[j];

  LOG_DBG3(lbm_pde) << "d_collision_step_d_f: index = " << index << " dfeqdux=" << dfeqdux.ToString() ;

  StdVector<double>& collMatrix = adjCollisions[index]; // storage that we hand over to LBM solver for adjoint simulation

  //partial derivative of collision operator with respect to f: CORRECT (FDM)
  for (unsigned int i = 0; i < n_q_; i++)
  {
    for (unsigned int j = 0; j < n_q_; j++)
    {
      block[i][j] = (i == j) * (1 - omega_) + omega_ * dfeqdf[i][j];
      collMatrix[GetMatrixElemId(j,i,n_q_)] = block[i][j]; // transposing included
    }
  }
}

void LatticeBoltzmannPDE::d_bounceback_d_f(int index, Matrix<double>& block)
{
  StdVector<double>& adjBB = adjCollisions[index];
  adjBB.Init();
  // bounce-back sensitivities; gradient is just a permutation matrix
  for (unsigned int dir = 0; dir < n_q_; dir++) {
    block[dir][invPDFDirections_[dir]] = 1.0;
    adjBB[GetMatrixElemId(invPDFDirections_[dir],dir,n_q_)] = 1.0;
  }
}

void LatticeBoltzmannPDE::d_inflow_d_f(int index, Matrix<double>& block, StdVector<double>& weight)
{
  // gradient of the velocity inlet boundary with respect to the design variables
  StdVector<double> dfeqdrho(n_q_);
  double dot, norm;

  if (dim_ == 2) assert(u_max_z_ == 0);

  LatticeBoltzmann::PDFDirectionVector transform;
  for (unsigned int i = 0; i < n_q_; i++)
  {
    transform = microVelDirections_[i];
    if (!parabolicInflow_) {
      dot = transform.off_x * u_max_x_ + transform.off_y * u_max_y_ + transform.off_z * u_max_z_;
      norm = u_max_x_ * u_max_x_ + u_max_y_ * u_max_y_ + u_max_z_ * u_max_z_;
    }
    else {
      // find index of vector u_in_ for this element
      int id = inlet.Find(index);
      dot = transform.off_x * u_in_[id][0] + transform.off_y * u_in_[id][1] + transform.off_z * u_in_[id][2];
      norm = u_in_[id][0] * u_in_[id][0] + u_in_[id][1] * u_in_[id][1] + u_in_[id][2] * u_in_[id][2]; // only one component is not zero
      LOG_DBG3(lbm_pde) << "d_inflow_d_f index=" << index << " i=" << i << " us1=" << u_in_[id][0] << " us2=" << u_in_[id][1] << " us3= " << u_in_[id][2] << " dot=" << dot << " norm=" << norm;
    }
    dfeqdrho[i] = weight[i] * (1. + 3 * dot + 4.5 * dot * dot - 1.5 * norm);
  }
   LOG_DBG3(lbm_pde) << "d_inflow_d_f index=" << index << " dfeqdrho=" << StdVector<double>::ToString(9, &dfeqdrho[0]);

  StdVector<double>& adjInlet = adjCollisions[index];

  for (unsigned int i = 0; i < n_q_; i++)
    for (unsigned int j = 0; j < n_q_; j++) {
      block[i][j] = dfeqdrho[i];
      adjInlet[GetMatrixElemId(j,i,n_q_)] = dfeqdrho[i];
    }
}

void LatticeBoltzmannPDE::d_outflow_d_f(int index, Matrix<double>& block, StdVector<double>& ux, StdVector<double>& uy, StdVector<double>& uz, StdVector<double>& dloc, StdVector<double>& weight)
{
  // gradient of the density outlet boundary condition with respect to the design variables
  StdVector<double> dfeqdux(n_q_), dfeqduy(n_q_), dfeqduz(n_q_);
  LatticeBoltzmann::PDFDirectionVector transform;
  double scale = 1.; // no design
  double density = 1.0;  // we have no other case
  double us1, us2, us3, dot;

  for (unsigned int i = 0; i < n_q_; i++)
  {
    transform = microVelDirections_[i];
    us1 = scale * ux[index];
    us2 = scale * uy[index];
    us3 = scale * uz[index];

    dot = transform.off_x * us1 + transform.off_y * us2 + transform.off_z * us3;
    dfeqdux[i] = weight[i] * density * (3 * transform.off_x + 9. * transform.off_x * dot - 3 * us1);
    dfeqduy[i] = weight[i] * density * (3 * transform.off_y + 9. * transform.off_y * dot - 3 * us2);
    dfeqduz[i] = weight[i] * density * (3 * transform.off_z + 9. * transform.off_z * dot - 3 * us3);
  }
  //gradient of u_x/u_y with respect to f
  StdVector<double> duxdf(n_q_), duydf(n_q_), duzdf(n_q_);

  double invdloc = 1.0 / dloc[index];

  for (unsigned int i = 0; i < n_q_ ; i++) {
    transform = microVelDirections_[i];
    duxdf[i] = invdloc * (-ux[index] + transform.off_x);
    duydf[i] = invdloc * (-uy[index] + transform.off_y);
    duzdf[i] = invdloc * (-uz[index] + transform.off_z);
  }

  StdVector<double>& adjOutlet = adjCollisions[index];
  //gradient of u_x with respect to f
  for (unsigned int i = 0; i < n_q_; i++)
    for (unsigned int j = 0;j < n_q_; j++) {
      block[i][j] = dfeqdux[i] * duxdf[j] + dfeqduy[i] * duydf[j] + dfeqduz[i] * duzdf[j];
      adjOutlet[GetMatrixElemId(j,i,n_q_)] = block[i][j];
    }

}

/*
 * 3 case have to be considered for each distribution function f_*:
 * 1. f_* corresponds to an element that is not on the boundary --> f_* influences only its neighbour
 * 2. f_* corresponds to a boundary element and is normal to the corresponding boundary edge
 *    --> functions on the boundary are not being collided, only propagated
 *    --> after propagation, the previous value of f_* now influences the new value of f_* at the same position and the new value of f_* at the neighbour position
 * 3. f_* corresponds to a boundary element and is tangent to the corresponding boundary edge --> equivalent to case 1
 */
void LatticeBoltzmannPDE::d_propagate_d_f(mapped_matrix<double>& Jprop, const mapped_matrix<double>& J)
{
  //gradient of the propagations step with respect to the design variables by Georg Pingen, University of Colorado (CU), Boulder, Colorado
  unsigned int x, y, z;
  int rows1,rows2;
  mapped_matrix<double>::const_iterator1 iter;
  LatticeBoltzmann::PDFDirectionVector transform;

  for(z = 0; z < n_z_ ; z++) {
    for(y = 0; y < n_y_ ; y++) {
      for(x = 0; x < n_x_ ; x++) {
        for (unsigned int dir = 0; dir < n_q_; dir++) {
          transform = microVelDirections_[dir];

          // distributions pointing outside the domain don't influence the simulation --> d_propagate/d_f = 0
          if (!PointsToBoundary(x,y,z,dir)) {
            rows1 = GetPdfIndex(x,y,z,dir);
            rows2 = GetPdfIndex(x + transform.off_x,y + transform.off_y, z + transform.off_z, dir);

            // case 1; if inverse distribution function of f_* points inside the domain, f_* is tangent to its boundary edge
            if (!PointsToBoundary(x,y,z,(invPDFDirections_)[dir])) {
              iter = J.find1(0, rows2, 0);
              for(mapped_matrix<double>::const_iterator2 it = iter.begin(); it != iter.end(); ++it) {
                Jprop(rows1,it.index2()) = J(rows2,it.index2());
              }
            }
            // case 2
            else {
              iter = J.find1(0, rows1, 0);
              for(mapped_matrix<double>::const_iterator2 it = iter.begin(); it != iter.end(); ++it) {
                Jprop(rows1,it.index2()) = J(rows1,it.index2());
              }
              iter = J.find1(0, rows2, 0);
              for(mapped_matrix<double>::const_iterator2 it = iter.begin(); it != iter.end(); ++it) {
                Jprop(rows1,it.index2()) += J(rows2,it.index2());
              }
            }
          }

        }
      }
    }
  }
}

Vector<double> LatticeBoltzmannPDE::d_pressuredrop_d_f(StdVector<double>& ux, StdVector<double>& uy, StdVector<double>& uz)
{
  // Calculation of gradient of pressure drop with respect to design variable; By Georg Pingen, University of Colorado (CU), Boulder, Colorado
//  StdVector<double> dUX(n_q_), dUY(n_q_), dUZ(n_q_);
  double dUX, dUY, dUZ;

  mapped_matrix<double> dPD(n_elems * n_q_, 1, n_elems * n_q_);
  double inletSize_inv = 1.0 / (double) inlet.GetSize();
  double outletSize_inv = 1.0 / (double) outlet.GetSize();
  double one_third = 1.0 / 3.0;

  // temporal storage for adjoint SRT
  StdVector<double> d_PD_d_f(n_q_);
  int index;

  for (unsigned int i = 0; i < inlet.GetSize(); i++) {
    index = inlet[i];
    for (unsigned int dir = 0; dir < n_q_; dir++) {
      dUX = microVelDirections_[dir].off_x - ux[index]; // inlet velocities are calculated from steady state solution and not the prescribed ones from xml file
      dUY = microVelDirections_[dir].off_y - uy[index];
      dUZ = microVelDirections_[dir].off_z - uz[index];
      dPD(GetPdfIndex(index,dir),0) = (one_third + 0.5 * (ux[index] * ux[index] + uy[index] * uy[index]) + uz[index] * uz[index] + ux[index] * dUX + uy[index] * dUY + uz[index] * dUZ) * inletSize_inv;
//      d_PD_d_f[dir] = dPD(GetPdfIndex(index,dir),0) * 1e4;
      d_PD_d_f[dir] = dPD(GetPdfIndex(index,dir),0);
    }
    d_pdrop_d_f[index] = d_PD_d_f;
  }

  for (unsigned int i = 0; i < outlet.GetSize(); i++) {
    index = outlet[i];
    for (unsigned int dir = 0; dir < n_q_; dir++) {
      dUX = microVelDirections_[dir].off_x - ux[index];
      dUY = microVelDirections_[dir].off_y - uy[index];
      dUZ = microVelDirections_[dir].off_z - uz[index];
      dPD(GetPdfIndex(index,dir),0) = -(one_third + 0.5 * (ux[index]*ux[index] + uy[index] * uy[index] + uz[index] * uz[index]) + ux[index] * dUX + uy[index] * dUY + uz[index] * dUZ) * outletSize_inv;
//      d_PD_d_f[dir] = dPD(GetPdfIndex(index,dir),0) * 1e4;
      d_PD_d_f[dir] = dPD(GetPdfIndex(index,dir),0);
    }
    d_pdrop_d_f[index] = d_PD_d_f;
  }

  Vector<double> rhs(n_elems * n_q_);

//  dPD *= 1e4;

  mapped_matrix<double> dFdf(n_elems * n_q_, 1, n_elems * n_q_);
  d_propagate_d_f(dFdf,dPD);
  for(unsigned int i = 0, n = rhs.GetSize(); i < n; i++)
    rhs[i] = dFdf(i,0);

  return rhs; // no copy constructor
}

Vector<Double> LatticeBoltzmannPDE::CalcVelocities(unsigned int idx)
{
  Vector<Double> velo;
  velo.Resize(dim_);
  double density = CalcLBMDensity(idx);
  velo[0] =  CalcVelocityX(idx, density);
  velo[1] = CalcVelocityY(idx, density);
  if (dim_ == 3)
    velo[2] = CalcVelocityZ(idx, density);
  return velo;
}

double LatticeBoltzmannPDE::CalcPressure(unsigned int idx) const
{
  double density = CalcLBMDensity(idx);
  double ux     = CalcVelocityX(idx, density);
  double uy     = CalcVelocityY(idx, density);
  double uz     = CalcVelocityZ(idx, density);

  if (dim_ == 2)
    assert(uz == 0);

  LOG_DBG2(lbm_pde) << "CP: " << "idx = " << idx << " density = " << density << " ux = " << ux << " uy = " << uy << " uz = " << uz << " p = " << density / 3.0 + 0.5 * density * (ux * ux + uy * uy + uz * uz);
  return density / 3.0 + 0.5 * density * (ux * ux + uy * uy + uz * uz);
}

double LatticeBoltzmannPDE::CalcPressureDrop()
{
  // Calculation of the pressure drop by Pingen
  double in = 0.0;
  for(unsigned int i = 0; i < inlet.GetSize(); i++) {
    in += CalcPressure(inlet[i]);
    LOG_DBG2(lbm_pde) << "CPD: idx=" << inlet[i] << " p=" << CalcPressure(inlet[i]);
  }

  double out = 0.0;
  for(unsigned int i = 0; i < outlet.GetSize(); i++) {
    out += CalcPressure(outlet[i]);
    LOG_DBG2(lbm_pde) << "CPD: idx=" <<outlet[i] << " p=" << CalcPressure(outlet[i]);
  }

  LOG_DBG2(lbm_pde) << "CPD: dP = " << in / inlet.GetSize() - out / outlet.GetSize();
//  return (in / inlet.GetSize() - out / outlet.GetSize()) * 1e4;
  return (in / inlet.GetSize() - out / outlet.GetSize());
}

Vector<Double> LatticeBoltzmannPDE::ExtractDistribution(unsigned int idx){
  Vector<double> solVec;
  solVec.Resize(n_q_);
  for (unsigned int i = 0; i < n_q_; i++) {
    solVec[i] = GetPdf(idx,i);
  }
  LOG_DBG3(lbm_pde) << "ED: " << " idx=" << idx << " pdfs=" << solVec.ToString();
  return solVec;
}

// ***********************************************************************
//   Obtain information on desired output quantities from parameter file
// ***********************************************************************
void LatticeBoltzmannPDE::ReadSpecialResults()
{
}

void LatticeBoltzmannPDE::DefinePrimaryResults() {
  //setting some dummy state variable
  shared_ptr<ResultInfo> res1( new ResultInfo);
  res1->resultType = TEST_DOF;

  res1->dofNames = "";
  res1->unit = "?";
  res1->definedOn = ResultInfo::MapSolTypeToDefinedOn(TEST_DOF);
  res1->entryType = ResultInfo::SCALAR;
  feFunctions_[TEST_DOF]->SetResultInfo(res1);
  results_.Push_back( res1 );
  availResults_.insert( res1 );
  res1->SetFeFunction(feFunctions_[TEST_DOF]);
  DefineFieldResult( feFunctions_[TEST_DOF], res1 );
}


void LatticeBoltzmannPDE::DefinePostProcResults() {

  shared_ptr<BaseFeFunction> feFct = feFunctions_[TEST_DOF];

  shared_ptr<ResultInfo> prob(new ResultInfo);
  prob->resultType = LBM_PROBABILITY_DISTRIBUTION;
  StdVector<std::string> probDofNames(n_q_); // "f_0", "f_1", ....
  for(unsigned int i=0, n = n_q_; i < n; i++ )
    probDofNames[i] = "f_" + boost::lexical_cast<std::string>(i);
  prob->dofNames = probDofNames;
  prob->unit = "";
  prob->entryType = ResultInfo::VECTOR;
  prob->definedOn = ResultInfo::MapSolTypeToDefinedOn(LBM_PROBABILITY_DISTRIBUTION);
  shared_ptr<CoefFunctionFormBased> probFunc;
  probFunc.reset(new CoefFunctionLBM<Double>(this,feFct, prob));
  DefineFieldResult(probFunc,prob);

  // =====================================================================
  // set solution information
  // =====================================================================
  PtrParamNode resultsNode = myParam_->Get("storeResults", ParamNode::PASS );

  // === solution by LBM solver ===
  // result LBM_PROBABILITY_DISTRIBUTION already defined in DefinePrimaryResults()
  shared_ptr<ResultInfo> dens(new ResultInfo);
  dens->resultType = LBM_DENSITY;
  dens->unit =  "";
  dens->dofNames = "";
  dens->entryType = ResultInfo::SCALAR;
  dens->definedOn = ResultInfo::MapSolTypeToDefinedOn(LBM_DENSITY);
  shared_ptr<CoefFunctionFormBased> densFunc;
  densFunc.reset(new CoefFunctionLBM<Double>(this,feFct, dens));
  DefineFieldResult(densFunc,dens);

  shared_ptr<ResultInfo> press(new ResultInfo);
  press->resultType = LBM_PRESSURE;
  press->unit =  "";
  press->dofNames = "";
  press->entryType = ResultInfo::SCALAR;
  press->definedOn = ResultInfo::MapSolTypeToDefinedOn(LBM_PRESSURE);
  shared_ptr<CoefFunctionFormBased> pressFunc;
  pressFunc.reset(new CoefFunctionLBM<Double>(this,feFct, press));
  DefineFieldResult(pressFunc,press);

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
  velo->definedOn = ResultInfo::MapSolTypeToDefinedOn(LBM_VELOCITY);
  shared_ptr<CoefFunctionFormBased> veloFunc;
  veloFunc.reset(new CoefFunctionLBM<Double>(this,feFct, velo));
  DefineFieldResult(veloFunc,velo);

  // === PSEUDO DENSITY for SIMP ===
  shared_ptr<ResultInfo> mechPD(new ResultInfo);
  mechPD->resultType = MECH_PSEUDO_DENSITY;
  mechPD->dofNames = "";
  mechPD->unit = "";
  mechPD->entryType = ResultInfo::SCALAR;
  mechPD->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_PSEUDO_DENSITY);
  mechPD->fromOptimization = true;
  shared_ptr<CoefFunctionFormBased> mechPDFunc;
  mechPDFunc.reset(new CoefFunctionLBM<Double>(this,feFct, mechPD));
  DefineFieldResult(mechPDFunc,mechPD);

  shared_ptr<ResultInfo> physicalPD(new ResultInfo);
  physicalPD->resultType = LBM_PHYSICAL_PSEUDO_DENSITY;
  physicalPD->dofNames = "";
  physicalPD->unit = "";
  physicalPD->entryType = ResultInfo::SCALAR;
  physicalPD->definedOn = ResultInfo::MapSolTypeToDefinedOn(LBM_PHYSICAL_PSEUDO_DENSITY);
  physicalPD->fromOptimization = true;
  shared_ptr<CoefFunctionFormBased> physicalPDFunc;
  physicalPDFunc.reset(new CoefFunctionLBM<Double>(this,feFct, physicalPD));
  DefineFieldResult(physicalPDFunc,physicalPD);
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
      if (n_q_ == 9) {
        if(!(ss >>  GetPdf(i,0) >>  GetPdf(i,1) >>  GetPdf(i,2) >>  GetPdf(i,3) >>  GetPdf(i,4) >>  GetPdf(i,5) >>  GetPdf(i,6) >>  GetPdf(i,7) >>  GetPdf(i,8)))
          EXCEPTION("error reading nine values in line " << (i+1) << " of file " << filename);
      }
      else if (n_q_ == 19){
        if(!(ss >>  GetPdf(i,0) >>  GetPdf(i,1) >>  GetPdf(i,2) >>  GetPdf(i,3) >>  GetPdf(i,4) >>  GetPdf(i,5) >>  GetPdf(i,6) >>  GetPdf(i,7) >>  GetPdf(i,8) >> GetPdf(i,9)
            >> GetPdf(i,10) >> GetPdf(i,11) >> GetPdf(i,12) >> GetPdf(i,13) >> GetPdf(i,14) >> GetPdf(i,15) >> GetPdf(i,16) >> GetPdf(i,17) >> GetPdf(i,18)))
          EXCEPTION("error reading nineteen values in line " << (i+1) << " of file " << filename);
      }
      i++;
    }
  }

  if(i != n_elems)
    EXCEPTION("read " << i << " lines in " << filename << " but expected " << n_elems);


  file.close();

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
    file <<  u_max_x_ << " " << u_max_y_ << " " << u_max_z_ << std::endl;

  file.close();

  std::cout << "++ CFS2LBM.dat created" << std::endl;
}

//might be necessary for debugging
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
      // out: -1 bb, -2 inlet, -3 outlet, 0 ... 1 porosity
      // out: 1 bb1, 2 inlet, 3 outlet, 0 for porosity
      int v = (int) -1.0 * elements[n_x_ * y + x];
      obst << std::max(0, v) << " "; // org porosity would be -1
      obst_m(y, x) = std::max(0, v);  // origin of the matrix is left upper!
    }
    obst << std::endl;
  }
  obst.close();

  LOG_DBG2(lbm_pde) << "EMF: obst matrix (!= obst.dat):\n" << obst_m.ToString();

  std::ofstream ns("non_sing.dat");
  // add one to be one based
  for(unsigned int i = 0, n = non_sing.GetSize(); i < n; i++)
    ns << (non_sing[i] + 1) << " " << std::endl; // space to allow diff with matlab's non_sing.dat
  ns.close();

  std::ofstream data("data.dat");
  data << n_x_ << std::endl;
  data << n_y_ << std::endl;
  data << 1.0 << std::endl; // penalty parameter
  data << u_max_x_ << std::endl;
  data << u_max_y_ << std::endl;
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

  // find out which elements are inlet, outlet, boundary and obstacle ones
  StdVector<Elem*> elems;
  StdVector<Elem*> boundaries;
  StdVector<Elem*> obst;
  // auxiliary vector
  // vector initialized with density value of inner cells
  elements.Resize(n_elems, 0.0);

  DesignSpace* space = domain->GetDesign();
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
        double val = space->GetErsatzMaterialFactor(idx, App::LBM);
        // the ordering of the design elements and the LBM elements might be different as it is defined for LBM
        // and might be arbitrary in the mesh which defines the ordering in the optimization design
        elements[elem_to_idx[elem->elemNum]] = val;
        rel.Push_back(elem_to_idx[elem->elemNum]); // storing id of nodes associated with design elements;
      }
    }
  }

  StdVector<std::string> regionNames;
  grd->GetRegionNames(regionNames);
  //specify regionId for boundary and inner region
  int boundId = -1;
  int innerId = -1;
  int inclusionId = -1;

  for (UInt i = 0; i<regionNames.GetSize(); ++i) {
    if (regionNames[i] == "boundary")
      boundId = i;
    else if (regionNames[i] == "obstacle")
      inclusionId = i;
    else
      innerId = i;
  }
  if (boundId == -1 || innerId == -1)
    EXCEPTION("Mesh must contain 2 regions: 'boundary' and another region");

  grd->GetElems(boundaries,boundId);
  grd->GetElems(elems,innerId);

  for (unsigned int i = 0; i < boundaries.GetSize(); ++i) {
    elements[boundaries[i]->elemNum - 1] = -1.0;
    if (!inlet.Contains(boundaries[i]->elemNum - 1) && !outlet.Contains(boundaries[i]->elemNum - 1))
      bb.Push_back(boundaries[i]->elemNum - 1); // if boundary node is neither inlet nor outlet node, then it's a bounce back node
  }
  // if there are obstacles in the domain
  if (inclusionId != -1) {
    grd->GetElems(obst,inclusionId);
    for (unsigned int i = 0; i < obst.GetSize(); ++i) {
      elements[obst[i]->elemNum - 1] = 1.0;
      obstacles.Push_back(obst[i]->elemNum - 1);
    }
  }
  for(unsigned int i = 0; i < inlet.GetSize(); ++i)
    elements[inlet[i]] = -2.0;
  for(unsigned int i = 0; i < outlet.GetSize(); ++i)
    elements[outlet[i]] = -3.0;
}

void LatticeBoltzmannPDE::SetupParabolicInflow()
{
  // for a parabolic profile only one prescribed velocity component is allowed to be inequal 0
  if (((u_max_x_ == 0.0) + (u_max_y_ == 0.0) + (u_max_z_ == 0.0)) != 2)
    EXCEPTION("For a parabolic flow profile only one inflow velocity component is allowed to be not 0!")
    u_in_.Resize(inlet.GetSize());

  int flow_index;
  double u_vertex;

  if (u_max_x_ != 0.0) {
    flow_index = 0;
    u_vertex = u_max_x_;
  }
  else if (u_max_y_ != 0.0) {
    flow_index = 1;
    u_vertex = u_max_y_;
  }
  else {
    flow_index = 2;
    u_vertex = u_max_z_;
  }
  // assume indices in inlet vector are ordered!
  int n_inlet = inlet.GetSize();

  unsigned int y_max = 0;
  unsigned int y_min = 9999;
  unsigned int z_max = 0;
  unsigned int z_min = 9999;

  StdVector< StdVector<unsigned int> > coordsInlet; // store coordinates of all inlet elements
  coordsInlet.Resize(n_inlet);
  for (int i = 0; i < n_inlet; ++i) {
    StdVector<unsigned int> tmp;
    tmp.Resize(3);
    tmp.Init(0);
    GetCoordinates(inlet[i],tmp[0],tmp[1],tmp[2]);
    if (tmp[1] > y_max) // find biggest inlet coordinate in y and z direction
      y_max = tmp[1];
    if (tmp[2] > z_max)
      z_max = tmp[2];
    if (tmp[1] < y_min) // find smallest inlet coordinate in y and z direction
      y_min = tmp[1];
    if (tmp[2] < z_min)
      z_min = tmp[2];
    coordsInlet[i] = tmp;
  }

  double vertex_y = y_min + (y_max - y_min) / 2.0;
  double vertex_z = 0;
  if (dim_ == 3)
    vertex_z = z_min + (z_max - z_min) / 2.0;
  double r_y = y_max - vertex_y + 1;
  double r_z = z_max - vertex_z + 1;

  for (int i = 0; i < n_inlet; ++i) {

    double disty = (double)coordsInlet[i][1] - vertex_y;

    double distz = (double)coordsInlet[i][2] - vertex_z;
    StdVector<double> tmp_u;
    tmp_u.Resize(3);
    tmp_u.Init(0.0);
    tmp_u[flow_index] = u_vertex * (1 - disty * disty / (r_y*r_y)) * (1 - distz * distz / (r_z*r_z));
    u_in_[i] = tmp_u;
  }
}

void LatticeBoltzmannPDE::SetNonSingualrityIndices()
{
  assert(non_sing.IsEmpty());
  assert(!elements.IsEmpty());

  non_sing.Reserve(n_q_ * n_x_ * n_y_ * n_z_);

  for(unsigned int z = 0; z < n_z_; z++)
  {
    for(unsigned int y = 0; y < n_y_; y++)
    {
      for(unsigned int x = 0; x < n_x_; x++)
      {
        unsigned int base = GetPdfIndex(x,y,z,Q_0); // 2D
        //  LOG_DBG3(lbm_pde) << "INSI: test y=" << y << " x=" << x << " k=" << k << " base=" << base << " val -> " << val;
        if(elements[GetIndex(x,y,z)] != LBM_NODE_TYPE_BB) // no bounce back
        {
          // fluid node distribution functions are not deleted
          for(unsigned int h = 0; h < n_q_; h++) // 2D
            non_sing.Push_back(base + h);
        }
        else
        {
          // outpointing boundary distributions are not inserted and
          // distributions which point towards a bounce-back boundary point

          for (unsigned int dir = 0; dir < n_q_; dir++) {
            LatticeBoltzmann::PDFDirectionVector f = microVelDirections_[dir];
            // test, if the node in direction f is not bounce back
            if (!PointsToBoundary(x,y,z,dir) && elements[GetIndex(x+f.off_x, y+f.off_y, z+f.off_z)] != LBM_NODE_TYPE_BB)
            {
              non_sing.Push_back(base + dir);
            }
          }
        }
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
  for (unsigned int i = 0; i < n_elems;i++)
  {
    for (unsigned int j = 0; j < n_q_;j++) {
      f << GetPdf(i,j)<<" ";
    }
    f << std::endl;
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
  for(mapped_matrix<double>::const_iterator1 it = M.begin1(); it != M.end1(); ++it)
  {
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
  f << "% Matrix extracted from boost::mapped_matrix\n";
  f << "%\n";
  f << M.size1() << "\t" << M.size2() << "\t" << number <<  "\n";
  f << ss.str() << std::endl;
  f.close();
}

// automatic testing for PointsToBoundary(...)
void LatticeBoltzmannPDE::TestPointsToBoundary()
{
  assert(!PointsToBoundary(0,0,0,Q_0));
  assert(!PointsToBoundary(0,0,0,Q_E));
  assert(!PointsToBoundary(0,0,0,Q_N));
  assert(PointsToBoundary(0,0,0,Q_W));
  assert(PointsToBoundary(0,0,0,Q_S));
  assert(!PointsToBoundary(0,0,0,Q_NE));
  assert(PointsToBoundary(0,0,0,Q_NW));
  assert(PointsToBoundary(0,0,0,Q_SW));
}

} // end of namespace CoupledField
