// =====================================================================================
// 
//       Filename:  buInt.cc
// 
//    Description:  Implementation file for BUIntegrators
//                  TAKE care:
//                  This file is just for inclusion in the header file!
// 
//        Version:  1.0
//        Created:  11/03/2011 08:08:21 PM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Andreas Hueppe (AHU), andreas.hueppe@uni-klu.ac.at
//        Company:  Universitaet Klagenfurt
// 
// =====================================================================================

namespace CoupledField
{

  template <class VEC_DATA_TYPE, bool SURFACE>
  BUIntegrator<VEC_DATA_TYPE, SURFACE>::
      BUIntegrator(BaseBOperator *bOp,
                   VEC_DATA_TYPE factor,
                   shared_ptr<CoefFunction> rhsCoef,
                   bool coordUpdate,
                   bool fullEvaluation,
                   bool extractReal,
                   const string &id)
      : LinearForm(coordUpdate),
        fullEvaluation_(fullEvaluation)
  {
    factor_ = factor;
    this->name_ = "RhsBUIntegrator";
    this->bOperator_ = bOp;

    this->rhsCoefs_ = rhsCoef;
    extractReal_ = extractReal;
    id_ = id;
    useVolume4Edge_ = false;
  }

  template <class VEC_DATA_TYPE, bool SURFACE>
  BUIntegrator<VEC_DATA_TYPE, SURFACE>::
      BUIntegrator(BaseBOperator *bOp,
                   VEC_DATA_TYPE factor,
                   shared_ptr<CoefFunction> rhsCoef,
                   const std::set<RegionIdType> &volRegions,
                   bool coordUpdate,
                   bool fullEvaluation,
                   bool extractReal,
                   const string &id,
                   bool useVolume4Edge)
      : LinearForm(coordUpdate),
        fullEvaluation_(fullEvaluation)
  {
    factor_ = factor;
    this->name_ = "RhsBUIntegrator";
    this->bOperator_ = bOp;
    extractReal_ = extractReal;
    id_ = id;
    useVolume4Edge_ = useVolume4Edge;

    assert(rhsCoef->GetDimType() == CoefFunction::VECTOR ||
           rhsCoef->GetDimType() == CoefFunction::SCALAR);
#ifndef NDEBUG
  if (rhsCoef->GetDimType() != CoefFunction::VECTOR &&
      rhsCoef->GetDimType() != CoefFunction::SCALAR)
  {
    Exception("BDB integrator expects the coefficient function to be vectorial or scalar!");
  }
#endif
  this->rhsCoefs_ = rhsCoef;
  volRegions_ = volRegions;
}

  template <class VEC_DATA_TYPE, bool SURFACE>
  void BUIntegrator<VEC_DATA_TYPE, SURFACE>::
      CalcElemVector(Vector<VEC_DATA_TYPE> &elemVec,
                     EntityIterator &ent)
  {

    assert(rhsCoefs_->GetDimType() != CoefFunction::NO_DIM);

    // Declare necessary variables
    const Elem *ptElem = ent.GetElem();
    Matrix<Double> bMat;
    Vector<VEC_DATA_TYPE> cVec;
    StdVector<LocPoint> intPoints;
    StdVector<Double> weights;
    UInt nrFncs = 0;
    VEC_DATA_TYPE fac(0.0);

    // Surface: inverse of jacobian
    Vector<VEC_DATA_TYPE> pt1;
    Vector<VEC_DATA_TYPE> pt2;
    Matrix<Double> JacT;
    Matrix<Double> TF;
    Matrix<Double> TFinv;

    // Obtain FE element from feSpace and integration scheme
    IntegOrder order;
    IntScheme::IntegMethod method;
    BaseFE *ptFe = ptFeSpace_->GetFe(ent, method, order);
    BaseFE *ptFeVol = NULL;
    if (useVolume4Edge_)
    {
      // call the specialized version when using edge elements and apply
      // the operator to the volume shape functions (altough the integral
      // is a surface integral)
      // But first we need to get the volume-element side that was specified in the
      // PDE definition (volRegions)
      RegionIdType volRegionId1;
      RegionIdType volRegionId2;
      if(ent.GetSurfElem()->ptVolElems.size() > 0){
        ent.GetGrid()->GetElemRegion(ent.GetSurfElem()->ptVolElems[0]->elemNum, volRegionId1);
      }
      if ((ent.GetSurfElem()->ptVolElems.size() > 1) && (ent.GetSurfElem()->ptVolElems[1] != NULL)){
        ent.GetGrid()->GetElemRegion(ent.GetSurfElem()->ptVolElems[1]->elemNum, volRegionId2);
      }

      bool isInVol1 = this->volRegions_.find(volRegionId1) != this->volRegions_.end();
      bool isInVol2 = this->volRegions_.find(volRegionId2) != this->volRegions_.end();

      if(isInVol1){
        this->surfaceVolReg_ = 0;
        ptFeVol = ptFeSpace_->GetFe(ent.GetSurfElem()->ptVolElems[0]->elemNum);
      }else if(isInVol2){
        this->surfaceVolReg_ = 1;
        ptFeVol = ptFeSpace_->GetFe(ent.GetSurfElem()->ptVolElems[1]->elemNum);
      }else{
        EXCEPTION("Specified surface region is not in the provided volume regions!");
      }
    }

    // nrFncs = ptFe->GetNumFncs();
  
    if (useVolume4Edge_)
    {
      nrFncs = ptFeVol->GetNumFncs();
    }
    else
    {
      nrFncs = ptFe->GetNumFncs();
    }  

    // Get shape map from grid
    shared_ptr<ElemShapeMap> esm =
        ent.GetGrid()->GetElemShapeMap(ptElem, this->coordUpdate_);

    // Get integration points
    intScheme_->GetIntPoints(Elem::GetShapeType(ptElem->type), method, order,
                             intPoints, weights);

    LocPointMapped lp;
    elemVec.Resize(nrFncs * Bdim_);
    elemVec.Init();

    // Pre-evaluate coefficient function in case of reduced accuracy
    if (!fullEvaluation_)
    {
      const ElemShape sh = Elem::shapes[ptElem->type];
      lp.Set(sh.midPointCoord, esm, sh.volume);
      if (rhsCoefs_->GetDimType() == CoefFunction::SCALAR)
      {
        cVec.Resize(1);
        rhsCoefs_->GetScalar(cVec[0], lp);
      }
      else
      {
        rhsCoefs_->GetVector(cVec, lp);
      }
    }

    // Loop over all integration points
    for (UInt i = 0; i < intPoints.GetSize(); i++)
    {

      // Calculate for each integration point the LocPointMapped
      if (SURFACE)
      {
        lp.Set(intPoints[i], esm, volRegions_, weights[i]);
      }
      else
      {
        lp.Set(intPoints[i], esm, weights[i]);
      }

      // calc factor
      fac = VEC_DATA_TYPE(lp.jacDet * weights[i]);
      fac *= factor_;

      if (useVolume4Edge_)
      {
        bOperator_->CalcOpMatTransposed(bMat, lp, ptFeVol);
      }
      else
      {
        bOperator_->CalcOpMatTransposed(bMat, lp, ptFe);
      }
      // Evaluate coefficient function in integration point
      // ( in case of full order)
      if (fullEvaluation_)
      {
        if (rhsCoefs_->GetDimType() == CoefFunction::SCALAR)
        {
          cVec.Resize(1);
          rhsCoefs_->GetScalar(cVec[0], lp);
        }
        else
        {
          rhsCoefs_->GetVector(cVec, lp);
          if (SURFACE && (ptFeSpace_->GetSpaceType() == FeSpace::HCURL))
          {
            // uxn
            pt1 = lp.normal;
            cVec.CrossProduct(pt1, pt2);

            // Jacobian of surface element
            lp.jac.Transpose(JacT);

            // Metric and its inverse
            TF = JacT * lp.jac;
            TF.Invert(TFinv);

            // Transformation of a function in curl space (see Zaglmayer Lemma 4.15)
            cVec = (TFinv * JacT) * pt2;
            fac *= lp.lpmVol->jacDet;
          }
        }
      }
      elemVec += bMat * cVec * fac;
    }
  }
}