// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_PREISACH_2004
#define FILE_PREISACH_2004

#include <list>

#include "MatVec/vector.hh"
#include "MatVec/matrix.hh"
#include "hysteresis.hh"

namespace CoupledField {

  class Preisach : public Hysteresis
  {
  public:
    //! constructor
    Preisach(Integer numElem, Double xSat, Double ysat, 
	     Matrix<Double>& preisachWeight, bool isVirgin);

    //!
    virtual ~Preisach();

    //!
    Double computeValue(Double& xVal, Integer idxElem, bool overwrite = true);

    //!
    Double computeInverseValue(Double xVal, Integer idxElem);

    //!
    Double computeValueAndUpdate(Double xVal, Integer idxElem, 
                                 bool overwrite = true);

    //! 
    Double getValue(  Integer idxElem);

    //! 
    UInt getStringLength( Integer idxElem ) {
      return  StringLenght_[ idxElem];
    };

    //!
    virtual Double GetIncX() {
      return dH_;
    };

    //!
    Double updateMinMaxList(Double newX, Integer idxElem, 
                          bool overwrite);

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
    Double normalizeInput(Double xInput);

    //!
    Double normalizeOutput(Double xInput);


    //! compute preisach weights;
    void computePreisachWeights();

  protected:

  private:

    Double Xsaturated_;
    Double YSaturated_;
    Double YRemnant_;

    bool isVirgin_;
    Integer actElem_;

    Vector<Double> lastVal_;
    Vector<Double> preisachSum_;

    Vector<Double>* strings_;
    Vector<Double>* helpStrings_;

    Vector<UInt> StringLenght_;
    UInt maxStringLength_;

    Matrix<Double> preisachWeights_;

    Double dH_; 
    Double eps_;
  };


} //end of namespace


#endif

