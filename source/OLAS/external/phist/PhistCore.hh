/*
 * PhistCore.hh
 *
 *  Created on: Mar 5, 2018
 *      Author: sri
 */

#ifndef OLAS_EXTERNAL_PHIST_PHISTCORE_HH_
#define OLAS_EXTERNAL_PHIST_PHISTCORE_HH_

#include <def_expl_templ_inst.hh>

#include <limits>
#include <string.h>

#include "MatVec/StdMatrix.hh"
#include "MatVec/CRS_Matrix.hh"
#include "MatVec/SCRS_Matrix.hh"
#include "MatVec/generatematvec.hh"
#include "MatVec/Matrix.hh"

#include "OLAS/algsys/SolStrategy.hh"
#include "OLAS/precond/generateprecond.hh"
#include "OLAS/precond/BasePrecond.hh"
#include "OLAS/solver/generatesolver.hh"
#include "OLAS/solver/BaseSolver.hh"
#include "OLAS/solver/BaseEigenSolver.hh"

#include "Utils/Timer.hh"
#include "Domain/Domain.hh"
#include "Driver/BaseDriver.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

#include "phist_enums.h"
#include "phist_config.h"
#include "phist_types.hpp"
#include "phist_jadaOpts.h"
#include "phist_jada.hpp"
#include "phist_void_aliases.h"
#include "phist_kernels.hpp"
#include "phist_core.hpp"
#include "phist_krylov.hpp"
#include "phist_precon.hpp"

namespace CoupledField {

class BaseMatrix;
class BaseVector;
class Flags;
class sparseMat_t; // phist matrix type
class StdMatrix;

class PhistCore{
public:
  PhistCore();
  sparseMat_t* InitMatrix(const BaseMatrix& cfs, sparseMat_t** phist, double scale);
  ~PhistCore();
  typedef struct {
    /** either stiff, mass, or damping */
    const StdMatrix* mat = NULL;
    /** scale value, e.g. to scale the B-Mat by 1/B[0,0]. Controlled by scale_mass*/
    double scale = 1.0;
  }SparseMatRowFuncService;


private:

  phist_comm_ptr comm_ = NULL;

};

} /* namespace CoupledField */

#endif /* OLAS_EXTERNAL_PHIST_PHISTCORE_HH_ */
