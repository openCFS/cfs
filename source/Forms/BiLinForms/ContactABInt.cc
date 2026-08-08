// ================================================================================================
/*!
 *       \file     ContactABInt.cc
 *       \brief    Penalty contact stiffness for one coupling block of a contact pair.
 *                 See ContactABInt.hh for the derivation and the block signs.
 *       \date     2026
 */
//================================================================================================

#include "ContactABInt.hh"

#include "Domain/ElemMapping/EntityLists.hh"
#include "Forms/ContactPointInterpolation.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

DEFINE_LOG(contactabint, "contactabint")

namespace CoupledField {

ContactABInt::ContactABInt(ContactInterface* iface,
                           BiLinearForm::CouplingDirection direction,
                           Mode mode)
  : BiLinearForm(),
    iface_(iface),
    direction_(direction),
    mode_(mode)
{
  this->name_ = (mode == MODE_STABILIZATION) ? "ContactStabABInt"
              : (mode == MODE_FRICTION)      ? "ContactFrictionABInt"
                                             : "ContactABInt";
}

ContactABInt::ContactABInt(const ContactABInt& right)
  : BiLinearForm(right),
    iface_(right.iface_),
    direction_(right.direction_),
    mode_(right.mode_)
{
}

void ContactABInt::SetFeSpace(shared_ptr<FeSpace> feSpace1, shared_ptr<FeSpace> feSpace2) {
  this->ptFeSpace1_ = feSpace1;
  this->ptFeSpace2_ = feSpace2;
}

// ================================================================================================

void ContactABInt::CalcElementMatrix(Matrix<Double>& elemMat,
                                     EntityIterator& ent1,
                                     EntityIterator& ent2) {

  const NcSurfElem* ncElem = ent1.GetNcSurfElem();
  const ContactNcSurfElem* ce = dynamic_cast<const ContactNcSurfElem*>(ncElem);
  if (ce == nullptr || ce->segment == nullptr) {
    EXCEPTION("ContactABInt expects a ContactNcSurfElem carrying a contact segment.");
  }

  const ContactSegment& seg = *(ce->segment);
  const UInt dim = iface_->GetGrid()->GetDim();

  // Which side supplies the test functions (rows) and which the trial functions (columns).
  // Must mirror SurfaceBiLinFormContext::MapEqns() for the same CouplingDirection:
  //   PRIM_* -> rows on the primary body, SEC_* -> rows on the secondary body.
  const bool rowIsPrimary = (direction_ == BiLinearForm::PRIM_PRIM ||
                             direction_ == BiLinearForm::PRIM_SEC);
  const bool colIsPrimary = (direction_ == BiLinearForm::PRIM_PRIM ||
                             direction_ == BiLinearForm::SEC_PRIM);

  // N_D = [ N_s , -N_m ], so a primary factor carries a minus sign.
  const Double rowSign = rowIsPrimary ? -1.0 : 1.0;
  const Double colSign = colIsPrimary ? -1.0 : 1.0;

  bool sized = false;

  // Only the points of THIS assembler element, i.e. those sharing its primary element.
  // Iterating the whole segment would assemble points belonging to a different primary
  // element into these equations -- see the header of ContactNcSurfElem.
  for (UInt k = 0; k < ce->pointIdx.GetSize(); ++k) {

    const ContactPoint& cp = seg.points[ce->pointIdx[k]];
    if (!cp.isProjected) {
      continue;
    }

    const Double gap = iface_->EvalCurrentGap(cp);
    const bool active = iface_->IsGapActive(seg, cp, gap);
    const bool wantActive = (mode_ != MODE_STABILIZATION);
    if (active != wantActive) {
      continue;
    }
    // g == 0 counts as active. It carries no force, but it must carry stiffness, otherwise
    // two bodies meshed exactly touching start from a singular system.

    // rows
    SurfElem* rowElem = rowIsPrimary ? cp.primElem : cp.secElem;
    const LocPoint& rowLp = rowIsPrimary ? cp.primLocal : cp.secLocal;
    CalcContactInterpolationMatrix(iface_->GetGrid(), rowElem, rowLp,
                                   rowIsPrimary ? iface_->GetPrimaryVolRegions()
                                                : iface_->GetSecondaryVolRegions(),
                                   rowIsPrimary ? this->ptFeSpace1_ : this->ptFeSpace2_,
                                   bMatRow_);

    // columns
    SurfElem* colElem = colIsPrimary ? cp.primElem : cp.secElem;
    const LocPoint& colLp = colIsPrimary ? cp.primLocal : cp.secLocal;
    CalcContactInterpolationMatrix(iface_->GetGrid(), colElem, colLp,
                                   colIsPrimary ? iface_->GetPrimaryVolRegions()
                                                : iface_->GetSecondaryVolRegions(),
                                   colIsPrimary ? this->ptFeSpace1_ : this->ptFeSpace2_,
                                   bMatCol_);

    const UInt nRow = bMatRow_.GetNumCols();
    const UInt nCol = bMatCol_.GetNumCols();

    if (!sized) {
      elemMat.Resize(nRow, nCol);
      elemMat.Init();
      sized = true;
    }

    if (mode_ == MODE_FRICTION) {

      // ------------------------------------------------------------------------------
      //  Consistent tangent of the Coulomb return map, both branches.
      //    K_RC = rowSign * colSign * INT w N_R^T D N_C dG ,   D = dt_T/d(u_s - u_m)
      // ------------------------------------------------------------------------------
      iface_->EvalCurrentSlip(cp, slip_);
      const bool stick = iface_->EvalFrictionTraction(seg, cp, gap, slip_, tracT_);

      // slip_ becomes the ELASTIC slip s = g_T - g_T^p, which is what the return map
      // linearizes about.
      for (UInt d = 0; d < dim; ++d) {
        slip_[d] -= cp.slipPlastic[d];
      }

      const Double epsT = iface_->GetTangentialPenalty(seg);
      const Double mu   = iface_->GetFrictionCoefficient();

      dMat_.Resize(dim, dim);
      dMat_.Init();

      if (stick) {
        // D = eps_T P_T, symmetric and positive semi-definite.
        for (UInt d = 0; d < dim; ++d) {
          for (UInt e = 0; e < dim; ++e) {
            dMat_[d][e] = epsT * (((d == e) ? 1.0 : 0.0) - cp.normal[d] * cp.normal[e]);
          }
        }
      } else {
        Double normS = 0.0;
        for (UInt d = 0; d < dim; ++d) {
          normS += slip_[d] * slip_[d];
        }
        normS = std::sqrt(normS);
        if (normS <= 0.0) {
          // Unreachable: the slip branch is entered only when ||eps_T s|| exceeds mu*p >= 0,
          // which requires ||s|| > 0. Skip rather than divide, so a future change to the
          // return map cannot turn this into a NaN.
          continue;
        }

        // The Coulomb limit itself depends on u through the normal pressure. Differentiate
        // the SAME expression the traction uses: t_N = min(lambda + eps_N g, 0), so the
        // derivative is -eps_N where the point actually carries pressure and 0 at the kink.
        const Double tN = cp.lambda + seg.penalty * gap;
        const Double dpdg =
            (tN < 0.0 && iface_->GetFrictionTangent()
                           == ContactInterface::FRICTION_TANGENT_CONSISTENT)
            ? -seg.penalty : 0.0;                        // dp/dg_N
        const Double p = -std::min(tN, 0.0);
        const Double radial = (normS > 0.0) ? (mu * p / normS) : 0.0;

        for (UInt d = 0; d < dim; ++d) {
          const Double ed = slip_[d] / normS;
          for (UInt e = 0; e < dim; ++e) {
            const Double ee = slip_[e] / normS;
            Double v = mu * dpdg * ed * cp.normal[e];
            v += radial * ((((d == e) ? 1.0 : 0.0) - cp.normal[d] * cp.normal[e]) - ed * ee);
            dMat_[d][e] = v;
          }
        }
      }

      const Double facF = rowSign * colSign * cp.weight;
      for (UInt i = 0; i < nRow; ++i) {
        for (UInt j = 0; j < nCol; ++j) {
          Double v = 0.0;
          for (UInt d = 0; d < dim; ++d) {
            for (UInt e = 0; e < dim; ++e) {
              v += bMatRow_[d][i] * dMat_[d][e] * bMatCol_[e][j];
            }
          }
          elemMat[i][j] += facF * v;
        }
      }
      continue;
    }

    // Project both interpolation matrices onto the contact normal:
    //   a_i = (N_row^T n)_i ,  b_j = (n^T N_col)_j
    nRow_.Resize(nRow);
    nRow_.Init();
    for (UInt i = 0; i < nRow; ++i) {
      Double v = 0.0;
      for (UInt d = 0; d < dim; ++d) {
        v += bMatRow_[d][i] * cp.normal[d];
      }
      nRow_[i] = v;
    }
    nCol_.Resize(nCol);
    nCol_.Init();
    for (UInt j = 0; j < nCol; ++j) {
      Double v = 0.0;
      for (UInt d = 0; d < dim; ++d) {
        v += bMatCol_[d][j] * cp.normal[d];
      }
      nCol_[j] = v;
    }

    const Double stiffness = (mode_ == MODE_CONTACT) ? seg.penalty
                                                     : iface_->GetStabilizationStiffness(seg);
    if (stiffness == 0.0) {
      continue;
    }
    const Double fac = rowSign * colSign * stiffness * cp.weight;

    for (UInt i = 0; i < nRow; ++i) {
      for (UInt j = 0; j < nCol; ++j) {
        elemMat[i][j] += fac * nRow_[i] * nCol_[j];
      }
    }
  }

  if (!sized) {
    // No active point on this segment. The assembler still needs a correctly shaped zero
    // matrix, so build it from the equation counts of both sides.
    StdVector<Integer> eqnRow, eqnCol;
    const Elem* rowVol = (rowIsPrimary ? ce->ptPrimary : ce->ptSecondary)->ptVolElems[0];
    const Elem* colVol = (colIsPrimary ? ce->ptPrimary : ce->ptSecondary)->ptVolElems[0];
    (rowIsPrimary ? this->ptFeSpace1_ : this->ptFeSpace2_)->GetElemEqns(eqnRow, rowVol);
    (colIsPrimary ? this->ptFeSpace1_ : this->ptFeSpace2_)->GetElemEqns(eqnCol, colVol);
    elemMat.Resize(eqnRow.GetSize(), eqnCol.GetSize());
    elemMat.Init();
  }

  LOG_DBG3(contactabint) << "dir=" << direction_ << " seg=" << ce->elemNum
      << " -> " << elemMat.GetNumRows() << "x" << elemMat.GetNumCols();
}

} // namespace CoupledField
