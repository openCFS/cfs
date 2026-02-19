#include "SGPHolder.hh"

#include <math.h>
#include <iterator>
#include <vector>
#include <algorithm>

#include "DataInOut/Logging/LogConfigurator.hh"
#include "Domain/CoefFunction/CoefFunctionOpt.hh"
#include "Driver/Assemble.hh"
#include "Forms/BiLinForms/BDBInt.hh"
#include "MatVec/Matrix.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Excitation.hh"
#include "Optimization/Optimization.hh"
#include "Optimization/OptimizationMaterial.hh"
#include "Optimization/TransferFunction.hh"
#include "Optimization/Optimizer/SGP.hh"
#include "PDE/StdPDE.hh"
#include "Utils/tools.hh"
#include "Utils/StdVector.hh"

DEFINE_LOG(sgpholder, "sgpholder")

using namespace CoupledField;
using std::pow;
using std::max;
using std::min;
using std::abs;
using std::string;
using boost::lexical_cast;

SGPHolder::SGPHolder(Optimization* opt, PtrParamNode pn) : BaseOptimizer(opt, pn, Optimization::SGP_SOLVER)
{
  space_ = optimization->GetDesign();
  n_elems_ = space_->GetNumberOfElements();

  int dim = domain->GetDim();
  sgp::InitParams initargs;
  initargs.dim = dim;

  Matrix<double> D = GetCoreStiffnessTensor();
  // copy to external sgp's data structure
  sgp::Matrix tens(sgp::n_VoigtNotation(dim));
  for (int i = 0; i < tens.rows(); i++)
    for (int j = 0; j < tens.cols(); j++)
      tens(i,j) = D[i][j];

  LOG_DBG3(sgpholder) << "stiffness of core material:" << D.ToString();

  initargs.D_core = tens;       // stiffness of core material

  // set FE volume for each element
  Vector<double> vol = GetElemVols();
  initargs.v_FE = sgp::Vector(vol.GetPointer(), vol.GetPointer() + vol.GetSize()); // element volumes

  sgp_ = new SGP(optimization, this, pn, initargs);
  opt_timer->Stop();
}

void SGPHolder::PostInitScale()
{
  BaseOptimizer::PostInitScale(1.0);
}
void SGPHolder::SolveProblem()
{
  assert(opt_timer->IsRunning());
  do {
    sgp_->Init();
    sgp_->SolveProblem();
  }
  while(restart_requested);
}

Matrix<double> SGPHolder::GetCoreStiffnessTensor()
{
  MaterialTensor<double> tens(VOIGT);
  DesignMaterial* dm = Optimization::context->dm;
  assert(dm != NULL);
  // in case we have multiple regions, but only one of them is subject to optimization
  // take first element in design domain
  StdVector<Elem*> elems;
  domain->GetGrid()->GetElems(elems,space_->GetRegionId());
  Elem* elem = elems[0];

  assert(elem->regionId == space_->GetRegionId());

  // FIXME: technical DesignElement::ALL_DESIGN is not correct, but it serves the cause
  // we want pure core stiffness tensor without density rho et.c -> provide external identity transfer function to design material
  dm->GetTensor(tens, DesignElement::ALL_DESIGNS, Optimization::context->stt, elem, DesignElement::NO_DERIVATIVE, VOIGT, true);

  return tens.GetMatrix(VOIGT);
}


Vector<double> SGPHolder::GetElemVols()
{
  Vector<double> vol(n_elems_);
  if(space_->IsRegular())
  { // we use 1/(the volume of the first element) as fraction
    double vole = domain->GetGrid()->GetElemShapeMap(space_->data[0].elem, false)->CalcVolume();
    vol.Init(vole);
  } else {
    for (unsigned int e = 0; e < n_elems_; e++) {
      DesignElement* de = &space_->data[e];
      vol[e] = de->CalcVolume();
    }
  } // if-else

  return vol;
}
