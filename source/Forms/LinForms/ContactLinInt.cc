// ================================================================================================
/*!
 *       \file     ContactLinInt.cc
 *       \brief    Reference-gap force of penalty contact.
 *                 See ContactLinInt.hh for the derivation and the side signs.
 *       \date     2026
 */
//================================================================================================

#include "ContactLinInt.hh"

#include "Domain/ElemMapping/EntityLists.hh"
#include "Forms/ContactPointInterpolation.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

DEFINE_LOG(contactlinint, "contactlinint")

namespace CoupledField {

ContactLinInt::ContactLinInt(ContactInterface* iface, bool rowIsPrimary, Mode mode)
  : LinearForm(),
    iface_(iface),
    rowIsPrimary_(rowIsPrimary),
    mode_(mode)
{
  this->name_ = (mode == MODE_FRICTION) ? "ContactFrictionLinInt" : "ContactLinInt";
}

ContactLinInt::ContactLinInt(const ContactLinInt& right)
  : LinearForm(right),
    iface_(right.iface_),
    rowIsPrimary_(right.rowIsPrimary_),
    mode_(right.mode_)
{
}

// ================================================================================================

void ContactLinInt::CalcElemVector(Vector<Double>& elemVec, EntityIterator& ent) {

  const NcSurfElem* ncElem = ent.GetNcSurfElem();
  const ContactNcSurfElem* ce = dynamic_cast<const ContactNcSurfElem*>(ncElem);
  if (ce == nullptr || ce->segment == nullptr) {
    EXCEPTION("ContactLinInt expects a ContactNcSurfElem carrying a contact segment.");
  }

  const ContactSegment& seg = *(ce->segment);
  const UInt dim = iface_->GetGrid()->GetDim();

  // N_D = [ N_s , -N_m ], so the primary side carries a minus sign. Identical to the rowSign
  // of ContactABInt, and it must stay identical: the two classes are the two halves of one
  // residual.
  const Double rowSign = rowIsPrimary_ ? -1.0 : 1.0;

  bool sized = false;

  for (UInt k = 0; k < ce->pointIdx.GetSize(); ++k) {

    const ContactPoint& cp = seg.points[ce->pointIdx[k]];
    if (!cp.isProjected) {
      continue;
    }

    const Double gap = iface_->EvalCurrentGap(cp);
    if (!iface_->IsGapActive(seg, cp, gap)) {
      continue;                       // separated: contributes nothing
    }

    SurfElem* rowElem = rowIsPrimary_ ? cp.primElem : cp.secElem;
    const LocPoint& rowLp = rowIsPrimary_ ? cp.primLocal : cp.secLocal;

    CalcContactInterpolationMatrix(iface_->GetGrid(), rowElem, rowLp,
                                   rowIsPrimary_ ? iface_->GetPrimaryVolRegions()
                                                 : iface_->GetSecondaryVolRegions(),
                                   this->ptFeSpace_,
                                   bMat_);

    const UInt nRow = bMat_.GetNumCols();

    if (!sized) {
      elemVec.Resize(nRow);
      elemVec.Init();
      sized = true;
    }

    if (mode_ == MODE_FRICTION) {

      // f_i += -rowSign * w * (N_row^T t_T)_i
      //
      // The COMPLETE tangential residual, not just a constant part of it: the friction law is
      // not quadratic in u, so unlike the normal term there is nothing a stiffness matrix
      // could manufacture. ContactABInt's MODE_FRICTION supplies the matching tangent as a
      // Newton form. See the header.
      iface_->EvalCurrentSlip(cp, slip_);
      iface_->EvalFrictionTraction(seg, cp, gap, slip_, tracT_);

      const Double facT = -rowSign * cp.weight;
      for (UInt i = 0; i < nRow; ++i) {
        Double v = 0.0;
        for (UInt d = 0; d < dim; ++d) {
          v += bMat_[d][i] * tracT_[d];
        }
        elemVec[i] += facT * v;
      }
      continue;
    }

    // f_i += -rowSign * w * (lambda + eps_N * g_0) * (N_row^T n)_i
    //
    // The augmented traction is  t_N = lambda + eps_N * (g_0 + N_D u)  on an active point, so
    // splitting off the part independent of u picks up lambda alongside eps_N * g_0. That is
    // the ONLY change augmented Lagrange makes to the assembled system: the stiffness K_c in
    // ContactABInt is untouched, because lambda is constant during a Newton solve.
    const Double fac = -rowSign * cp.weight * (cp.lambda + seg.penalty * cp.gap0);
    if (fac == 0.0) {
      continue;
    }

    for (UInt i = 0; i < nRow; ++i) {
      Double v = 0.0;
      for (UInt d = 0; d < dim; ++d) {
        v += bMat_[d][i] * cp.normal[d];
      }
      elemVec[i] += fac * v;
    }
  }

  if (!sized) {
    StdVector<Integer> eqn;
    const Elem* rowVol = (rowIsPrimary_ ? ce->ptPrimary : ce->ptSecondary)->ptVolElems[0];
    this->ptFeSpace_->GetElemEqns(eqn, rowVol);
    elemVec.Resize(eqn.GetSize());
    elemVec.Init();
  }

  LOG_DBG3(contactlinint) << "side=" << (rowIsPrimary_ ? "primary" : "secondary")
      << " seg=" << ce->elemNum << " -> " << elemVec.GetSize() << " entries";
}

} // namespace CoupledField
