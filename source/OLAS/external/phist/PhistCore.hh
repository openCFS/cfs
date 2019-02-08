/*
 * PhistCore.hh
 *
 *  Created on: Mar 5, 2018
 *      Author: sri
 *      Description: The Phist Core Class serves as intereface containing common functions of both the linear solver
 *       and eigen solver. The solvers share the matcreate functions and a few common operations.
 *       The linear and eigen solvers are derived from this and their respective base sovlers.
 *       This is done to avoid code duplications.
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

#include "Utils/Timer.hh"
#include "Domain/Domain.hh"
#include "Driver/BaseDriver.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

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
#include "phist_schur_decomp.h"


namespace CoupledField {

class BaseMatrix;
class BaseVector;
class Flags;
class sparseMat_t; // phist matrix type
class StdMatrix;

class PhistCore
{
public:
  PhistCore();
  ~PhistCore();

protected:

  /** creates the matrix. If the space is already occupied (!= NULL) it is deleted first
   * @param see it as referenece to the pointer we store the stuff to. It shall point to NULL in the beginning!
   * @return new content of phist */
  template<class TYPE>
  sparseMat_t* InitMatrix(const BaseMatrix& cfs, sparseMat_t** phist, double scale);

  template<class TYPE>
  static int SparseMatRowFunc(ghost_gidx row, ghost_lidx* row_nnz, ghost_gidx* row_col, void* values, void* service_void);

  typedef struct
  {
    /** either stiff, mass, or damping */
    const StdMatrix* mat = NULL;
    /** scale value, e.g. to scale the B-Mat by 1/B[0,0]. Controlled by scale_mass*/
    double scale = 1.0;
  } SparseMatRowFuncService;

private:

  // we dont't use the MPI features of phist and set it all to a minimum
  phist_comm_ptr comm_ = NULL;
};

} /* namespace CoupledField */

#endif /* OLAS_EXTERNAL_PHIST_PHISTCORE_HH_ */
