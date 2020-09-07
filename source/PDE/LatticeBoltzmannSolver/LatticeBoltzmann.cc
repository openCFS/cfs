#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
# ifdef _OPENMP
  #include <omp.h>
#endif

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ProgramOptions.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ResultHandler.hh"
#include "Domain/Domain.hh"
#include "Driver/BaseDriver.hh"
#include "Driver/StaticDriver.hh"
#include "Driver/Assemble.hh"
#include "Driver/SolveSteps/BaseSolveStep.hh"
#include "Driver/FormsContexts.hh"
#include "Optimization/Optimization.hh"
#include "PDE/LatticeBoltzmannSolver/LatticeBoltzmann.hh"
#include "Utils/Timer.hh"
#include "PDE/BasePDE.hh"


namespace CoupledField
{

using std::fstream;
using std::ios;

DEFINE_LOG(lbm, "lbm")

// instantiation of the static elements
Enum<LatticeBoltzmann::Direction>         LatticeBoltzmann::directions;

LatticeBoltzmann::LatticeBoltzmann(int dim, int sizeX, int sizeY, int sizeZ, double ux, double uy, double uz, StdVector< StdVector<double> > uin, double omega, int maxIterations, double maxTolerance, bool plot, int writeFrequency)
{
  assert(dim == 2 || dim == 3);
  // n_q_: number of discrete directions in this model, e.g. n_q_ for D3Qn_q_
  if (dim == 2) {
    assert(sizeZ == 1);
    n_q_ = 9;
  }
  else
    n_q_ = 19;

  dim_ = dim;
  sizeX_ = sizeX;
  sizeY_ = sizeY;
  sizeZ_ = sizeZ;
  ux_ = ux;
  uy_ = uy;
  uz_ = uz;
  omega_nu_ = omega; // this relaxation rate is directly related to the fluid's viscosity
  maxIter_ = maxIterations;
  maxTol_ = maxTolerance;
  writeFrequency_ = writeFrequency;
  numWriteResults_ = 0;

  nNodes_ = sizeX_ * sizeY_ * sizeZ_;

  u_in.Resize(nNodes_);
  u_in = uin;

  plot_ = plot;

  //matrix of the probability distributions
  LOG_DBG(lbm) << "Allocating arrays for " << nNodes_ * n_q_ << " PDFs (" << (sizeof(double) * nNodes_ * n_q_ * 2.0 / 1024.0 / 1024.0) << " MiB)";

  pdfs_.Resize(2);
  pdfs_[0].Resize(nNodes_ * n_q_);
  pdfs_[1].Resize(nNodes_ * n_q_);

  adjPdfs_.Resize(2);
  adjPdfs_[0].Resize(nNodes_ * n_q_);
  adjPdfs_[1].Resize(nNodes_ * n_q_);
  adjPdfs_[0].Init(0.0);
  adjPdfs_[1].Init(0.0);

  // microVelDirections stores information about the n_q_ microscopic velocities/directions of D3Qn_q_ model
  microVelDirections.Resize(n_q_);
  weights.Resize(n_q_);

  SetEnums();
  InitWeights();
  SetInvDirections();
  TestInvDirections();

  scales.Resize(nNodes_);
  rel.Resize(nNodes_);
  bb.Resize(2 * sizeX_ + 2 * sizeY_ + 2 * sizeZ_);

  cur_  = 0;
  next_ = 1;
  adjCur_ = 0;
  adjNext_ = 1;

  lbmCalls_ = 0;
  lbmAdjCalls_ = 0;

  SetMicroVelocities();

  //initlialize function pointers in dependence on problem's dimension
  if (dim == 2) {
    prop_coll_step = &LatticeBoltzmann::Prop_coll_step2D;
  }
  else
  {
    prop_coll_step = &LatticeBoltzmann::Prop_coll_step3D;
  }
  adjCollTimer_.reset(new Timer("adjColl",true));
  adjPropTimer_.reset(new Timer("adjProp",true));

  PtrParamNode LBMInfo = domain->GetInfoRoot()->Get(ParamNode::SUMMARY)->Get("LBMSolver", ParamNode::INSERT);
  LBMInfo->Get("timer", ParamNode::APPEND)->SetValue(adjCollTimer_);
  LBMInfo->Get("timer", ParamNode::APPEND)->SetValue(adjPropTimer_);
}

LatticeBoltzmann::~LatticeBoltzmann()
{
}

void LatticeBoltzmann::CalcVelocities(const Vector<double>& pdfs, double& ux, double& uy, double& uz)
{

  double density = CalcDensity(pdfs);
  double jx = 0;
  double jy = 0;
  double jz = 0;

  if (n_q_ == 9) {
    jx = (pdfs[Q_E] - pdfs[Q_W]) + (pdfs[Q_NE] - pdfs[Q_SW]) + (pdfs[Q_SE] - pdfs[Q_NW]);
    jy = (pdfs[Q_N] - pdfs[Q_S]) + (pdfs[Q_NE] - pdfs[Q_SW]) + (pdfs[Q_NW] - pdfs[Q_SE]);
  }
  else {
    jx = (pdfs[Q_E] - pdfs[Q_W]) + (pdfs[Q_NE] - pdfs[Q_SW]) + (pdfs[Q_SE] - pdfs[Q_NW]) + (pdfs[Q_TE] - pdfs[Q_BW]) + (pdfs[Q_BE] - pdfs[Q_TW]);
    jy = (pdfs[Q_N] - pdfs[Q_S]) + (pdfs[Q_NW] - pdfs[Q_SE]) + (pdfs[Q_NE] - pdfs[Q_SW]) + (pdfs[Q_TN] - pdfs[Q_BS]) + (pdfs[Q_BN] - pdfs[Q_TS]);
    jz = (pdfs[Q_T] - pdfs[Q_B]) + (pdfs[Q_TW] - pdfs[Q_BE]) + (pdfs[Q_TE] - pdfs[Q_BW]) + (pdfs[Q_TN] - pdfs[Q_BS]) + (pdfs[Q_TS] - pdfs[Q_BN]);
  }

  ux = jx / density;
  uy = jy / density;
  uz = jz / density;
}

StdVector<double>* LatticeBoltzmann::Iterate(const StdVector<double>& elements, PtrParamNode in)
{
  int it = 0;
  // counts number of written steps when not in optimization
  int count = 0;

  // this flag cannot be set in constructor, since information about optimization is not available when this constructor is called
  if (domain->GetOptimization() != NULL && domain->GetOptimization()->GetOptimizerType() != Optimization::EVALUATE_INITIAL_DESIGN) {
    writeIntermediateResults_ = false;
  }
  else
    writeIntermediateResults_ = true;

  double R = 1.0;
  bool steady_state = false;

  std::ofstream plot;
  if(plot_)
    plot.open(std::string(progOpts->GetSimName() + ".lbm.dat").c_str());

  InitializePdfs();
  SetupDataStructures(elements);

  assert((int) elements.GetSize() == nNodes_);
  #pragma ivdep
  for (int i = 0; i < nNodes_; ++i) {
    scales[i] = 1.0 - elements[i];

//    LOG_DBG3(lbm) << "Element " << i << " has density " << elements[i] << " and porosity " << scales[i];
  }

  Timer timer;
  timer.Start();

  in->Get("converged")->SetValue("running");

  while(it < maxIter_ && !steady_state && R <= 1000)
  {
    // -- Combined propagation and collision step -------------------------
    (this->*prop_coll_step)(cur_, next_);
    // -- Bounce back step ------------------------------------------------
    Prop_coll_bounce_back(next_); //same code for both 2D and 3D
    // -- Inlet condition -------------------------------------------------
    Prop_coll_velinlet(next_);
    // -- Outlet condition ------------------------------------------------
    Prop_coll_densoutlet(next_);

    if((it == 0 || it % 100 == 0)) // check convergence
    {
      R = CalcResidual(cur_,next_);

      if(R <= maxTol_)
        steady_state = true;

      if(plot_)
        plot << it << "\t" << R << "\n";

      in->Get("iterations")->SetValue(it);
      in->Get("residuum")->SetValue(R);
      domain->GetInfoRoot()->ToFile(); // is not written when called too often
    }

    cur_  = (cur_  + 1) % 2;
    next_ = (next_ + 1) % 2;

    it++;

    if (writeIntermediateResults_) {
      if (it % writeFrequency_ == 0) {
        domain->GetDriver()->StoreResults(count,(double) it);
        count++;
      }
    }
  }

  timer.Stop();

  PtrParamNode node = in->Get(ParamNode::PROCESS)->Get("forward", ParamNode::APPEND); // write out how many lbm iterations until convergence
  node->Get("number")->SetValue(lbmCalls_);
  node->Get("iterations")->SetValue(it);
  node->Get("residuum")->SetValue(R);
  node->Get("converged")->SetValue(steady_state);

  if(plot_) {
    plot << it << "\t" << R << "\n";
    plot.flush();
  }

  double wt = timer.GetWallTime();
  double performance;
  if (dim_ == 3)
    performance = (sizeX_ - 1) * (sizeY_ - 1) * (sizeZ_ - 1) * it / wt / 1e6;
  else
    performance = (sizeX_ - 1) * (sizeY_ - 1) * it / wt / 1e6;

//  std::cout << "performance: " << (sizeX_ - 1) << "*" << (sizeY_ - 1) << "*" << (sizeZ_ - 1) << "*" << it << "/" << wt << "/1e6" << "=" << performance << std::endl;

  node->Get("wall")->SetValue(wt);
  node->Get("cpu")->SetValue(timer.GetCPUTime());
  node->Get("MFLUP_s")->SetValue(performance);
  in->Get("memory")->SetValue(sizeof(double) * nNodes_ * n_q_ * 2.0 / 1024.0 / 1024.0);

  if(R >= 1000)
    EXCEPTION("In LBM iteration " << it << " residuum " << R << " too large ... abort");

  if(!steady_state)
    EXCEPTION("internal LBM simulation could not converge: iterations: " << it << " residuum: " << R);

  numIterations_ = it;
  numWriteResults_ = count;

  lbmCalls_++; // first solver call is call number 0 (to match iteration numbering of optimizer)

  return &(pdfs_[cur_]);
}

void LatticeBoltzmann::InitializePdfs()
{
#pragma omp parallel for
  for (int elem = 0; elem < nNodes_; elem++) {
    #pragma ivdep
    for (int  dir = 0; dir < n_q_; dir++) {
      PDF(0, elem, dir) = weights[dir];
      PDF(1, elem, dir) = weights[dir];
    }
  }
}

void LatticeBoltzmann::InitializeAdjPdfs()
{
  adjPdfs_[0].Init(0.0);
  adjPdfs_[1].Init(0.0);
}

void LatticeBoltzmann::SetMicroVelocities()
{
  microVelDirections[Q_0] = PDFDirectionVector(0,0,0);
  microVelDirections[Q_E] = PDFDirectionVector(1,0,0);
  microVelDirections[Q_W] = PDFDirectionVector(-1,0,0);
  microVelDirections[Q_N] = PDFDirectionVector(0,1,0);
  microVelDirections[Q_S] = PDFDirectionVector(0,-1,0);
  microVelDirections[Q_NE] = PDFDirectionVector(1,1,0);
  microVelDirections[Q_NW] = PDFDirectionVector(-1,1,0);
  microVelDirections[Q_SE] = PDFDirectionVector(1,-1,0);
  microVelDirections[Q_SW] = PDFDirectionVector(-1,-1,0);
  if (dim_ == 3) {
    microVelDirections[Q_T] = PDFDirectionVector(0,0,1);
    microVelDirections[Q_B] = PDFDirectionVector(0,0,-1);
    microVelDirections[Q_TN] = PDFDirectionVector(0,1,1);
    microVelDirections[Q_TS] = PDFDirectionVector(0,-1,1);
    microVelDirections[Q_TE] = PDFDirectionVector(1,0,1);
    microVelDirections[Q_TW] = PDFDirectionVector(-1,0,1);
    microVelDirections[Q_BN] = PDFDirectionVector(0,1,-1);
    microVelDirections[Q_BS] = PDFDirectionVector(0,-1,-1);
    microVelDirections[Q_BE] = PDFDirectionVector(1,0,-1);
    microVelDirections[Q_BW] = PDFDirectionVector(-1,0,-1);
  }
}

void LatticeBoltzmann::SetEnums()
{
  directions.SetName("Qn_q_ directions");
  directions.Add(Q_0,"C");
  directions.Add(Q_E,"E");
  directions.Add(Q_N,"N");
  directions.Add(Q_W,"W");
  directions.Add(Q_S,"S");
  directions.Add(Q_T,"T");
  directions.Add(Q_B,"B");
  directions.Add(Q_NE,"NE");
  directions.Add(Q_NW,"NW");
  directions.Add(Q_SW,"SW");
  directions.Add(Q_SE,"SE");
  directions.Add(Q_TN,"TN");
  directions.Add(Q_BS,"BS");
  directions.Add(Q_TS,"TS");
  directions.Add(Q_BN,"BN");
  directions.Add(Q_TE,"TE");
  directions.Add(Q_BW,"BW");
  directions.Add(Q_TW,"TW");
  directions.Add(Q_BE,"BE");
}

void LatticeBoltzmann::InitWeights()
{
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

void LatticeBoltzmann::SetInvDirections()
{
  invPDFDirections.Resize(n_q_);

  invPDFDirections[Q_0] = Q_0;
  invPDFDirections[Q_E] = Q_W;
  invPDFDirections[Q_N] = Q_S;
  invPDFDirections[Q_W] = Q_E;
  invPDFDirections[Q_S] = Q_N;
  invPDFDirections[Q_NE] = Q_SW;
  invPDFDirections[Q_SW] = Q_NE;
  invPDFDirections[Q_NW] = Q_SE;
  invPDFDirections[Q_SE] = Q_NW;

  if (dim_ == 3) {
    invPDFDirections[Q_T] = Q_B;
    invPDFDirections[Q_B] = Q_T;
    invPDFDirections[Q_BS] = Q_TN;
    invPDFDirections[Q_TS] = Q_BN;
    invPDFDirections[Q_TN] = Q_BS;
    invPDFDirections[Q_BN] = Q_TS;
    invPDFDirections[Q_TE] = Q_BW;
    invPDFDirections[Q_BW] = Q_TE;
    invPDFDirections[Q_TW] = Q_BE;
    invPDFDirections[Q_BE] = Q_TW;
  }
}

void LatticeBoltzmann::TestInvDirections()
{
  assert(GetInvDirection(Q_0) == Q_0);
  assert(GetInvDirection(Q_E) == Q_W);
  assert(GetInvDirection(Q_N) == Q_S);
  assert(GetInvDirection(Q_W) == Q_E);
  assert(GetInvDirection(Q_S) == Q_N);
  assert(GetInvDirection(Q_NE) == Q_SW);
  assert(GetInvDirection(Q_NW) == Q_SE);
  assert(GetInvDirection(Q_SW) == Q_NE);
  assert(GetInvDirection(Q_SE) == Q_NW);

  if (dim_ == 3) {
    assert(GetInvDirection(Q_T) == Q_B);
    assert(GetInvDirection(Q_B) == Q_T);
    assert(GetInvDirection(Q_TN) == Q_BS);
    assert(GetInvDirection(Q_TS) == Q_BN);
    assert(GetInvDirection(Q_BS) == Q_TN);
    assert(GetInvDirection(Q_BN) == Q_TS);
    assert(GetInvDirection(Q_TE) == Q_BW);
    assert(GetInvDirection(Q_BW) == Q_TE);
    assert(GetInvDirection(Q_TW) == Q_BE);
    assert(GetInvDirection(Q_BE) == Q_TW);
  }
}

void LatticeBoltzmann::SetupDataStructures(const StdVector<double>& elements)
{
  rel.Clear();
  bb.Clear();
  inlet.Clear();
  outlet.Clear();

  for(int elem = 0; elem < nNodes_; elem++)
  {
    double porosity = elements[elem];

    if (LbmNodeTypeIsFluid(porosity))
      rel.Push_back(elem);
    else if (LbmNodeTypeIsBB(porosity))
      bb.Push_back(elem);
    else if (LbmNodeTypeIsInlet(porosity))
      inlet.Push_back(elem);
    else if (LbmNodeTypeIsOutlet(porosity))
      outlet.Push_back(elem);
    else if (LbmNodeIsObstacle(porosity))
      obst.Push_back(elem);
  }
}

std::string LatticeBoltzmann::ToString(const StdVector<StdVector<int> >& data)
{
  std::stringstream ss;

  for(unsigned int  e = 0; e < data.GetSize(); e++)
  {
    ss << e << ": (";
    for(unsigned int  d = 0; d < data[e].GetSize(); d++)
    {
      ss << data[e][d];
      if(d < data[e].GetSize() - 1)
        ss << ", ";
    }
    ss << ") ";
  }
  return ss.str();
}

void LatticeBoltzmann::CreateOutput(const char * file, int cur)
{
  // for debug purposes
  fstream f;
  f.precision(16);
  f.open(file, ios::out);

  for(int i = 0; i < nNodes_; i++) {
    for(int  j = 0; j < n_q_; j++) {
      f << PDF(cur, i, j) << " ";
    }

    f << std::endl;
  }
  f.close();
}


void LatticeBoltzmann::Prop_step()
{
  int tmp_x, tmp_y, tmp_z;
  // perform a propagation step
  for (int z = 0; z < sizeZ_; ++z) {
    for (int y = 0; y < sizeY_ ; ++y) {
      for (int x = 0; x < sizeX_ ; ++x) {
        #pragma ivdep
        for (int  dir = 0; dir < n_q_; dir++) {
          int invDir = GetInvDirection((Direction)dir);
          if (PointsToBoundary(x,y,z,invDir)) { // if the neighbor element that I want to access is outside the domain, keep current value
            tmp_x = x; // here we only set the coordinates
            tmp_y = y;
            tmp_z = z;
          }
          // else: standard propagation (get value from neighbor pdf)
          else {
            tmp_x = microVelDirections[invDir].off_x + x;
            tmp_y = microVelDirections[invDir].off_y + y;
            tmp_z = microVelDirections[invDir].off_z + z;
          }

          PDF(next_,x,y,z,dir) = PDF(cur_, tmp_x, tmp_y,  tmp_z, dir);
        }
      }
    }
  }

  pdfs_[cur_] = pdfs_[next_];
}

// Calculates macroscopic density for given pdfs of an element
double LatticeBoltzmann::CalcDensity(const Vector<double>& pdfs)
{
  double sum = 0;
  for (int  dir = 0; dir < n_q_; dir++) {
    sum += pdfs[dir];
  }
  return sum;
}

/************************************************** 2D operators *****************************************************/

void LatticeBoltzmann::Prop_coll_step2D(int cur, int next)
{
  int z = 0;
  #pragma omp parallel default(none) shared(next, cur, z)
  {
    Vector<double> pdfs;
    pdfs.Resize(n_q_);
    #pragma omp parallel for collapse(2)
    for (int y = 0; y < sizeY_ ; ++y) {
      for (int x = 0; x < sizeX_ ; ++x) {

        int index= GetIndex(x,y,z);

        // sum: macroscopic density is sum over all discrete distributions of an element
        double sum = 0;
        double tmp_ux = 0;
        double tmp_uy = 0;
        double tmp_us = 0;
        int tmp_x,tmp_y;
        double tmp;

        // propagation
        for (int  dir = 0; dir < n_q_; dir++) {
          int invDir = GetInvDirection((Direction)dir);
          if (PointsToBoundary(x,y,z,invDir)) { // if the neighbor element that I want to access is outside the domain, keep current value
            tmp_x = x; // here we only set the coordinates
            tmp_y = y;
          }
          // else: standard propagation (get value from neighbor pdf)
          else {
            tmp_x = microVelDirections[invDir].off_x + x;
            tmp_y = microVelDirections[invDir].off_y + y;
          }

          //store current pdf values in array for better accessing
          pdfs[dir] = PDF(cur, tmp_x, tmp_y,  z, dir); // accessed pdf value can be the old one or the neighbor's value
          sum += pdfs[dir];
          tmp_ux += microVelDirections[dir].off_x*pdfs[dir];
          tmp_uy += microVelDirections[dir].off_y*pdfs[dir];
        }

        // macroscopic scaling by design variable
        double scale = scales[index];

        tmp_ux = scale * tmp_ux / sum;
        tmp_uy = scale * tmp_uy / sum;
        tmp_us = 1.5 * (tmp_ux * tmp_ux + tmp_uy * tmp_uy);
        tmp_ux = 3.0 * tmp_ux;
        tmp_uy = 3.0 * tmp_uy;

        // propagation and collision in one step
        #pragma ivdep
        for (int  dir = 0; dir < n_q_; dir++)
        {
          tmp = microVelDirections[dir].off_x * tmp_ux + microVelDirections[dir].off_y * tmp_uy ;
          // no collision on the boundaries
          if (x == 0 || y == 0 || x == sizeX_ - 1 || y == sizeY_ - 1)
            PDF(next, x, y, z, dir) = pdfs[dir];
          else {
            PDF(next, x, y, z, dir) = pdfs[dir] + omega_nu_ * (sum * weights[dir]  * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us) - pdfs[dir]);
          }
        }

      }
    }
  }

}

void LatticeBoltzmann::Prop_coll_velinlet(int cur)
{
  #pragma omp parallel default(none) shared(cur)
  {
    Vector<double> pdfs(n_q_);
    #pragma omp for
    for(unsigned int  i = 0; i < inlet.GetSize(); i++) {
      int index = inlet[i];

      for (int  dir = 0; dir < n_q_; dir++) {
        pdfs[dir] = PDF(cur,index,dir);
      }

      double sum = CalcDensity(pdfs);

      double tmp_ux = ux_; // velocity at inlet is prescribed
      double tmp_uy = uy_;
      double tmp_uz = uz_;

      if (u_in.GetSize() != 0) {
        tmp_ux = u_in[i][0];
        tmp_uy = u_in[i][1];
        if (dim_ == 3)
          tmp_uz = u_in[i][2];
      }

      double tmp_us = 1.5 * (tmp_ux * tmp_ux + tmp_uy * tmp_uy + tmp_uz * tmp_uz);
      tmp_ux = 3.0 * tmp_ux;
      tmp_uy = 3.0 * tmp_uy;
      tmp_uz = 3.0 * tmp_uz;
      for (int  dir = 0; dir < n_q_; dir++) {
        double tmp = microVelDirections[dir].off_x * tmp_ux + microVelDirections[dir].off_y * tmp_uy + microVelDirections[dir].off_z * tmp_uz;
        PDF(cur, index, dir)  = sum * weights[dir]  * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us);
      }
    }
  }
}

//
// Performs a bounce back step.
//
void LatticeBoltzmann::Prop_coll_bounce_back(int cur)
{
  #pragma omp parallel default(none) shared(cur)
  {
    Vector<double> pdfs(n_q_);
    #pragma omp for
    for(unsigned int  i = 0; i < bb.GetSize(); i++) {
      int index = bb[i];

      for (int  dir = 0; dir < n_q_; dir++) {
        pdfs[dir] = PDF(cur, index, dir);
      }
      for (int  dir = 0; dir < n_q_; dir++) {
        PDF(cur, index, GetInvDirection((Direction)dir)) = pdfs[dir];
      }
    }
  }
}

//
// Density outlet condition.
//
void LatticeBoltzmann::Prop_coll_densoutlet(int cur)
{
  #pragma omp parallel default(none) shared(cur)
  {
    Vector<double> pdfs;
    pdfs.Resize(n_q_);

    #pragma omp for
    for(unsigned int  i = 0; i < outlet.GetSize(); i++) {
      int index = outlet[i];

      for (int  dir = 0; dir < n_q_; dir++) {
        pdfs[dir] = PDF(cur,index,dir);
      }

      double tmp_ux, tmp_uy, tmp_uz, tmp_us, sum, tmp;

      CalcVelocities(pdfs, tmp_ux, tmp_uy, tmp_uz);

      sum = 1.0; // the enforced density

      tmp_us = 1.5 * (tmp_ux * tmp_ux + tmp_uy * tmp_uy + tmp_uz * tmp_uz);
      tmp_ux = 3.0 * tmp_ux;
      tmp_uy = 3.0 * tmp_uy;
      tmp_uz = 3.0 * tmp_uz;
      for (int  dir = 0; dir < n_q_; dir++) {
        tmp = microVelDirections[dir].off_x * tmp_ux + microVelDirections[dir].off_y * tmp_uy + microVelDirections[dir].off_z * tmp_uz;
        PDF(cur, index, dir) =  sum * weights[dir] * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us);
      }
    }
  }
}


/************************************************** adjoint 2D operators *****************************************************/
void LatticeBoltzmann::AdjointPropagation(int next)
{
  PDFDirectionVector transform;
  #pragma omp parallel for default(none) private(transform) shared(next) collapse(3)
  for (int z = 0; z < sizeZ_; z++)
    for (int y = 0; y < sizeY_; y++)
      for (int x = 0; x < sizeX_; x++)
      {
        // adjoint propagation
        for (int  dir = 0; dir < n_q_; dir++) {
          transform = microVelDirections[dir];
          double value = 0.0;
          int id1 =  GetIndex(x,y,z);
          int rows1 = GetPdfIndex(id1,dir);
          int id2 = GetIndex(x + transform.off_x,y + transform.off_y, z + transform.off_z);
          int rows2 = GetPdfIndex(id2, dir);
          // distributions pointing outside the domain don't influence the simulation --> d_propagate/d_f = 0
          if (!PointsToBoundary(x,y,z,dir)) {
            // case 1: f_* corresponds to an element that is not on the boundary --> f_* influences only its neighbour
            if (!PointsToBoundary(x,y,z,(invPDFDirections)[dir])) {
              value = tmpPdfs_[rows2];
            }
            else { // case 2
              value = tmpPdfs_[rows1] + tmpPdfs_[rows2];
//              value += tmpPdfs_[rows2];
            }
          }
          APDF(next,id1,dir) = value;
        }
      }
}

StdVector<double>* LatticeBoltzmann::IterateAdjointSRT(PtrParamNode info,const StdVector<StdVector<double> >& collisionMatrices, const StdVector<StdVector<double> >& d_pdrop_d_f)
{
  tmpPdfs_.Resize(nNodes_ * n_q_);

  int it = 0;
  double R = 0.0;
  bool steady_state = false;

  InitializeAdjPdfs();

  std::ofstream plot;
  if(plot_)
    plot.open(std::string(progOpts->GetSimName() + ".adjLbm.dat").c_str());

  Timer timer;
  timer.Start();
  StdVector<double> pdfs, tmp;
  while(it < maxIter_ && !steady_state && R <= 1000)
  {
    adjCollTimer_->Start();
    #pragma omp parallel default(none) private(pdfs,tmp) shared(collisionMatrices,d_pdrop_d_f)
    {
      pdfs.Resize(n_q_);
      tmp.Resize(n_q_);
      /***************** Adjoint SRT collision ***/
      #pragma omp for schedule(static)
      for (int index = 0; index < nNodes_; index++)
      {
        #pragma ivdep
        for (int dir = 0; dir < n_q_; dir++)
          pdfs[dir] = APDF(adjCur_,index,dir);

        // adjoint collision: f* = d_pdrop_d_f + (d_coll_d_f)^T * f
//        Matrix<double> d_coll_d_f = collisionMatrices[index]; // collision matrices are already transposed
//        d_pd_d_f = d_pdrop_d_f[index];
//        d_coll_d_f.Mult(pdfs,tmp);
        MultLinMatrixVector(collisionMatrices[index],pdfs,tmp);

        #pragma ivdep
        for (int dir = 0; dir < n_q_; dir++) {
          tmpPdfs_[GetPdfIndex(index,dir)] = -d_pdrop_d_f[index][dir] + tmp[dir];
        }
      }
    }

    adjCollTimer_->Stop();

    adjPropTimer_->Start();
    AdjointPropagation(adjNext_);
    adjPropTimer_->Stop();

    if((it == 0 || it % 100 == 0))
    {
      R = CalcAdjResidual(adjCur_,adjNext_);
      if(R <= maxTol_)
        steady_state = true;
      if(plot_) {
        plot << it << "\t" << R << "\n";
        plot.flush();
      }
    }

    adjCur_  = (adjCur_  + 1) % 2;
    adjNext_ = (adjNext_ + 1) % 2;

    it++;
  }

  timer.Stop();

  PtrParamNode node = info->Get(ParamNode::PROCESS)->Get("adjoint", ParamNode::APPEND); // write out how many lbm iterations until convergence
  node->Get("number")->SetValue(lbmAdjCalls_);
  node->Get("iterations")->SetValue(it);
  node->Get("residuum")->SetValue(R);
  node->Get("converged")->SetValue(steady_state);

  double wt = timer.GetWallTime();
  double performance;

  if (dim_ == 3)
    performance = (sizeX_ - 1) * (sizeY_ - 1) * (sizeZ_ - 1) * it / wt / 1e6;
  else
    performance = (sizeX_ - 1) * (sizeY_ - 1) * it / wt / 1e6;

  node->Get("wall")->SetValue(wt);
  node->Get("cpu")->SetValue(timer.GetCPUTime());
  node->Get("MFLUP_s")->SetValue(performance);
  info->Get("memory")->SetValue(sizeof(double) * nNodes_ * n_q_ * 2.0 / 1024.0 / 1024.0);

  lbmAdjCalls_++;

  if(R >= 1000)
    EXCEPTION("In adjoint SRT iteration " << it << " residuum " << R << " too large ... abort");

  if(!steady_state)
    EXCEPTION("Adjoint SRT simulation could not converge: iterations: " << it << " residuum: " << R);

  return &(adjPdfs_[adjCur_]);
}

/************************************************** 3D operators *****************************************************/
void LatticeBoltzmann::Prop_coll_step3D(int cur, int next)
{

  int x, y, z = 0;
  double tmp, tmp_ux, tmp_uy, tmp_uz, tmp_us, scale, sum;

  int index;

  int tmp_x, tmp_y, tmp_z;
  StdVector<double> pdfs;

  #pragma omp parallel default(none)\
      private(index), \
      private(tmp_ux, tmp_uy, tmp_uz, tmp_us, scale, sum, tmp, x, y, z, tmp_x, tmp_y, tmp_z, pdfs), \
      shared(next, cur)
  {
    pdfs.Resize(n_q_);
    #pragma omp for schedule(static)//collapse(3)
    for (z = 0; z < sizeZ_; ++z) {
      for (y = 0; y < sizeY_; ++y) {
        for (x = 0; x < sizeX_; ++x) {
          index= GetIndex(x,y,z);

          // sum: macroscopic density is sum over all discrete distributions of an element
          sum = 0;
          tmp_ux = 0;
          tmp_uy = 0;
          tmp_uz = 0;
          tmp_us = 0;

          for (int  dir = 0; dir < n_q_; dir++) {
            int invDir = GetInvDirection((Direction)dir);
            if (PointsToBoundary(x,y,z,invDir)) { // boundary case
              tmp_x = x;
              tmp_y = y;
              tmp_z = z;
            }
            else { // standard propagation rule
              tmp_x = microVelDirections[invDir].off_x + x;
              tmp_y = microVelDirections[invDir].off_y + y;
              tmp_z = microVelDirections[invDir].off_z + z;
            }

            //store current pdf values in array for better accessing
            pdfs[dir] = PDF(cur, tmp_x, tmp_y, tmp_z, dir);
            sum += pdfs[dir];
            tmp_ux += microVelDirections[dir].off_x*pdfs[dir];
            tmp_uy += microVelDirections[dir].off_y*pdfs[dir];
            tmp_uz += microVelDirections[dir].off_z*pdfs[dir];
          }

          // macroscopic scaling by design variable
          scale = scales[index];

          tmp_ux = scale * tmp_ux / sum;
          tmp_uy = scale * tmp_uy / sum;
          tmp_uz = scale * tmp_uz / sum;
          tmp_us = 1.5 * (tmp_ux * tmp_ux + tmp_uy * tmp_uy + tmp_uz * tmp_uz);

          tmp_ux = 3.0 * tmp_ux;
          tmp_uy = 3.0 * tmp_uy;
          tmp_uz = 3.0 * tmp_uz;

          for (int  dir = 0; dir < n_q_; dir++) {
            tmp = microVelDirections[dir].off_x * tmp_ux + microVelDirections[dir].off_y * tmp_uy + microVelDirections[dir].off_z * tmp_uz;
            if (x == 0 || y == 0 || x == sizeX_ - 1 || y == sizeY_ - 1 || z == 0 || z == sizeZ_ - 1)
              PDF(next, index, dir) = pdfs[dir];
            else {
              PDF(next, index, dir) = pdfs[dir] + omega_nu_ * ( sum * weights[dir]  * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us) - pdfs[dir]);
            }
          }
        }
      }
    }
  }
}

} // end namespace
