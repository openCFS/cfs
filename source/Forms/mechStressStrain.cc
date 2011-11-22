// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include "Forms/mechStressStrain.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "PDE/mechPDE.hh"
#include "PDE/elecPDE.hh"
#include "MatVec/vector.hh"


DECLARE_LOG(ss)
DEFINE_LOG(ss, "stress_strain")

namespace CoupledField
{

  // =============================================================================
  // base class for mechanical stress computation
  // =============================================================================

  // base class for calculation of mechanical stresses
  template <class TYPE>
  MechStressStrain<TYPE>::MechStressStrain(BaseMaterial* matData, SubTensorType type)
    : linElastInt(matData, type)
  {
    name_ = "MechStressStrain";
  }


  template <class TYPE>
  MechStressStrain<TYPE>::~MechStressStrain()
  {
  }

  template <class TYPE>
  void MechStressStrain<TYPE>::CalcStressStrainResult(MechPDE* mech, shared_ptr<BaseResult> res, SolutionType st, bool density)
  {
    assert(st == MECH_STRESS || st == MECH_STRAIN || st == VON_MISES_STRAIN || st == VON_MISES_STRESS);

    bool do_von_mises = st == VON_MISES_STRAIN || st == VON_MISES_STRESS;

    LOG_DBG3(ss) << "MSS::CSSR: st=" << SolutionTypeEnum.ToString(st) << " d=" << density << " vM=" << do_von_mises;

    unsigned int sdim = mech->GetStressStrainDim();

    Vector<TYPE> vec;
    vec.Init(0);

    Result<TYPE> &  actRes = dynamic_cast<Result<TYPE>&>(*res);
    EntityIterator it = actRes.GetEntityList()->GetIterator();

    Vector<TYPE> & actVal = actRes.GetVector();
    actVal.Resize( actRes.GetEntityList()->GetSize() * (do_von_mises ? 1 : sdim) );

    // element solution
    Matrix<TYPE> sol;

    double jac_det = 0.0;
    double factor = 0.0;

    // for the density we need the coordinates
    Matrix<double> coords;
    // tmp and M for von Mises
    Vector<TYPE> tmp;
    const Matrix<double>& m = mech->GetVonMisesMatrix(sdim == 6 ? 3 : 2);

    // Fetch material: As we assume, that all elements belong to
    // one and the same region, we simply take the subdomain of the first
    // element
    it.Begin();

    SetAnsatzFct( mech->GetResultInfos()[0]->fctType );
    // loop over elements
    for ( it.Begin(); !it.IsEnd(); it++ )
    {
      // we integrate over the element by averages summation and then multiplying with the volume
      const Vector<Double>& intWeights = it.GetElem()->ptElem->GetIntWeights();

      // the real element volume
      mech->getPDE_grid()->GetElemNodesCoord(coords, it.GetElem()->connect, true); // updated coordinates
      double elem_vol = it.GetElem()->ptElem->CalcVolume(coords, isaxi_);

      // reset target such that we can sum up
      if(!do_von_mises)
        for(UInt iDof = 0; iDof < sdim; iDof++ )
          actVal[it.GetPos()*sdim + iDof] = 0.0;
      // to sum up von Mises inner product
      TYPE inner = 0.0;

      Vector<Double>* intPoints = it.GetElem()->ptElem->GetIntPoints(); // an array of Vectors :(
      LOG_DBG3(ss) << "MSS::CSSR: el=" << it.GetElem()->elemNum << " #ip=" << intPoints->GetSize()
                        << " NIP=" << it.GetElem()->ptElem->GetNumIntPoints() << " #w=" << intWeights.GetSize();

      assert(intWeights.GetSize() == it.GetElem()->ptElem->GetNumIntPoints());

      // loop over the integration points.
      Matrix<Double> elemCoord;
      mech->getPDE_grid()->GetElemNodesCoord( elemCoord, it.GetElem()->connect );

      //set element solution once
      mech->getPDESolution()->GetElemSolutionAsMatrix(sol, it);
      SetActElemSol(sol);
      for(UInt ip = 1, n = it.GetElem()->ptElem->GetNumIntPoints(); ip <= n; ip++) // is NOT intPoints->GetSize()!
      {
        SetIntPoint(intPoints[ip-1]); // fuck 1-based!!

        if(st == MECH_STRAIN || st == VON_MISES_STRAIN)
          CalcStrainVec(vec,ip,it);
        else
          CalcStressVec(vec,ip,it);

        // we need the Jacobi determinant for a good averaging.
        jac_det = it.GetElem()->ptElem->CalcJacobianDetAtIp(ip, coords, it.GetElem());

        // if we are axisymmetric, we have to account for the radial weighting,
        // i.e. we have to multiply the jacobian determinant by 2*pi*r(global x-coordinate of the point)
        // we need to compensate for the summation and handle the density (which is not! normalizing by element volume).
        if( isaxi_ ) {
          Vector<Double> globIp;
          it.GetElem()->ptElem->Local2GlobalCoord( globIp, intPoints[ip-1], elemCoord, it.GetElem());
          jac_det *= globIp[0]*2*PI;
        }

        factor = intWeights[ip-1] * jac_det / (density ? 1.0 : elem_vol); // fuck 1-based!
        UnsetIntPoint();

        LOG_DBG3(ss) << "MSS::CSSR: el=" << it.GetElem()->elemNum << " ip=" << ip << " vec=" << vec.ToString() << " jac_det=" << jac_det
            << " w=" << intWeights[ip-1] << " ev=" << elem_vol << " d=" << density << " f=" << factor;

        // this is the tensor case, we copy the summed up elements
        if(!do_von_mises)
        {
          for( UInt iDof = 0; iDof < sdim; iDof++ )
            actVal[it.GetPos()*sdim + iDof] += vec[iDof] * factor;
        }
        // this is the von Mises case. We sum up the inner product
        else
        {
          tmp = m * vec;
          TYPE squared = vec.Inner(tmp);
          inner += squared * factor; // factor needs to be outside the inner product!
          LOG_DBG3(ss) << "MSS::CSSR: vm: squared=" << squared << " inner=" << inner;
        }
      }

      // after integrating we take the norm of the von Mises result
      if(do_von_mises)
        actVal[it.GetPos()] = std::sqrt(inner);

      if(do_von_mises) {
        LOG_DBG3(ss) << "MSS::CSSR: vM " << " el=" << it.GetElem()->elemNum << " -> " << actVal[it.GetPos()];
      }
      else {
        LOG_DBG3(ss) << "MSS::CSSR: el=" << it.GetElem()->elemNum << " -> " << StdVector<TYPE>::ToString(sdim, actVal.GetPointer() + it.GetPos()*sdim);
      }
    }
  }



  /// calculates of stresses T (vector notation)
  // T = c . S with c the material tensor snd S the linear strains
  // S = B_lin * u
  // see Habil. M. Kaltenbacher

  template <class TYPE>
  void MechStressStrain<TYPE>::CalcStressVec(Vector<TYPE>& stressVec, UInt ip, EntityIterator& ent)
  {
    Matrix<TYPE> dMat;
    linElastInt::calcDMat(dMat, ent.GetElem()); // this includes the pseudo density for topology optimization!

    Vector<TYPE> linStrain;
    CalcStrainVec(linStrain, ip, ent);
    stressVec.Resize(linStrain.size());
    stressVec = dMat * linStrain;

    // LOG_DBG3(ss) << "MSS::CalcStressVec e=" << ent.GetElem()->elemNum << " ip=" << ip << " dMat=" << dMat.ToString();
    // LOG_DBG3(ss) << "MSS::CalcStressVec e=" << ent.GetElem()->elemNum << " ip=" << ip << " strain=" << linStrain.ToString() << " stress="
    //             << stressVec.ToString() << " u=" << displVec.ToString() << " B=" << linBMat.ToString(); // << " coord=" << ptCoord_.ToString();
  }

  /// calculates green-lagrangian strains (linear part, vector notation)
  // see Habil. M. Kaltenbacher

  template <class TYPE>
  void MechStressStrain<TYPE>::CalcStrainVec(Vector<TYPE>& strainVec, UInt ip, EntityIterator& ent)
  {
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );

    // convert displacement of all elem nodes into one vector:
    // (uNode1X, uNode1Y, uNode2X, uNode2Y, ...)
    Vector<TYPE> displVec;
    elemDisp_.ConvertToVec_AppendCols(displVec);

    // linear differential operator B_lin
    Matrix<Double> linBMat;
    linElastInt::CalcBMat(linBMat, ip, ptCoord_);
    strainVec.Resize(linBMat.GetNumRows());
    strainVec = linBMat * displVec;
  }



#ifdef __GNUC__
  // Explicite template instantiation
  template class MechStressStrain<Double>;
  template class MechStressStrain<Complex>;
#endif

} // end namespace CoupledField
