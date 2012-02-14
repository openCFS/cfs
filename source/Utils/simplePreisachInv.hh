// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_SIMPLPREISACHINV_HH
#define FILE_SIMPLPREISACHINV_HH


#include "General/defs.hh"
#include "MatVec/matrix.hh"
#include "MatVec/vector.hh"
#include "hysteresis.hh"

namespace CoupledField {

  class SimplePreisachInv : public Hysteresis
  {
  public:
    //! constructor
    SimplePreisachInv(Integer numElem, Double xSat, Double ysat, 
                      Matrix<Double>& preisachWeight, bool isVirgin);

    //!
    virtual ~SimplePreisachInv();

    //!
    Double computeValue(Double& yVal, Integer idxElem, 
                        bool overwrite = true);

    void SetPreviousYval( Double yval, Integer idxElem ) {
      previousYval_[idxElem] = normalizeYval( yval );
    }

    Double getActXval ( UInt idxElem ) {
      return actXval_[idxElem]*Xsaturated_;
    };

    //!
    void SetTimeStepVal(Double dt) 
    {;};

    //!
    Double EvalEverett(Double x1, Double x2, Integer idx);

    //!
    Double everett(Double x1, Double x2);

    //!
    Double everettPixel(Double x1, Double x2);

    //!
    Double normalizeXval(Double xInput);

    //!
    Double normalizeYval(Double& yInput);


    //! compute preisach weights;
    void computePreisachWeights();

  protected:

  private:

    Double Xsaturated_;
    Double Ysaturated_;
    Double YRemnant_;

    bool isVirgin_;
    Integer actElem_;

    Vector<Double> previousYval_;
    Vector<Double> actXval_;
    Vector<Double> preisachSum_;

    Vector<Double>* strings_;
    Vector<Double>* helpStrings_;

    StdVector<UInt> StringLenght_;
    UInt maxStringLength_;

    Matrix<Double> preisachWeights_;
    Vector<Double> epts_;
    UInt M_;
    Double eps_;
    Double nu0_;
  };


} //end of namespace


#endif

