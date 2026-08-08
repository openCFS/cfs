// ================================================================================================
/*!
 *       \file     ContactPointInterpolation.hh
 *       \brief    Interpolation operator of the volume element behind a contact point.
 *
 *                 Shared by the two contact integrators so that they cannot drift apart:
 *
 *                   ContactABInt   (Forms/BiLinForms)  -- the penalty stiffness  K_c
 *                   ContactLinInt  (Forms/LinForms)    -- the reference-gap force -eps*g_0*...
 *
 *                 Both halves of the same residual  r_c = -eps INT N_D^T n g_0 dG - K_c u  must
 *                 use the *same* shape functions at the *same* points, otherwise the constant
 *                 force and the tangent describe slightly different discrete gaps and Newton
 *                 stops converging quadratically. Keeping this in one inline function is the
 *                 cheapest way to guarantee that.
 *
 *       \date     2026
 */
//================================================================================================

#ifndef FILE_CFS_CONTACTPOINTINTERPOLATION_HH
#define FILE_CFS_CONTACTPOINTINTERPOLATION_HH

#include <set>

#include "General/Environment.hh"
#include "MatVec/Matrix.hh"
#include "Domain/ElemMapping/Elem.hh"
#include "Domain/ElemMapping/SurfElem.hh"
#include "Domain/ElemMapping/ElemShapeMap.hh"
#include "Domain/Mesh/Grid.hh"
#include "FeBasis/FeSpace.hh"
#include "FeBasis/BaseFE.hh"
#include "FeBasis/H1/H1Elems.hh"
#include "Forms/Operators/IdentityOperator.hh"

namespace CoupledField {

//! Identity (interpolation) operator of the volume element behind a surface point.

//! Both contact contexts (SurfaceBiLinFormContext, SurfaceLinFormContext) hand out the
//! equation numbers of the VOLUME element adjacent to the interface surface, not of the
//! surface element. The operator matrix must therefore be built on the volume element too,
//! evaluated at the volume-local image of the surface point -- which is what
//! LocPointMapped::SetWithSurface() produces in lpmVol. Same arrangement as
//! SurfaceNitscheABInt.
//!
inline void CalcContactInterpolationMatrix(Grid* grid,
                                           SurfElem* surfElem,
                                           const LocPoint& lp,
                                           const std::set<RegionIdType>& volRegs,
                                           shared_ptr<FeSpace> space,
                                           Matrix<Double>& bMat) {

  const UInt dim = grid->GetDim();

  shared_ptr<ElemShapeMap> esm = grid->GetElemShapeMap(surfElem, false);
  LocPointMapped lpm;
  lpm.SetWithSurface(lp, esm, volRegs, 0.0);

  const Elem* volElem = surfElem->ptVolElems[0];
  BaseFE* ptFe = space->GetFe(volElem->elemNum);

  if (dim == 3) {
    IdentityOperator<FeH1, 3, 3> op;
    op.CalcOpMat(bMat, *lpm.lpmVol, ptFe);
  }
  else {
    IdentityOperator<FeH1, 2, 2> op;
    op.CalcOpMat(bMat, *lpm.lpmVol, ptFe);
  }
}

} // namespace CoupledField

#endif /* FILE_CFS_CONTACTPOINTINTERPOLATION_HH */
