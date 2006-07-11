#ifndef FILE_PREISACH_2004
#define FILE_PREISACH_2004

#include <list>

#include "Utils/vector.hh"
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
    Double computeValue(Double xVal, Integer idxElem);

    //!
    void updateMinMaxList(Double newX, Integer idxElem);

    //!
    void SetTimeStepVal(Double dt) 
    {;};

    //!
    Double everett(Double x1, Double x2);

    //!
    Double everettPixel(Double x1, Double x2);

    //!
    Double normalizeInput(Double xInput);

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

    std::vector<Double>* strings_;
    std::vector<Double>* helpStrings_;

    std::vector<UInt> StringLenght_;
    UInt maxStringLength_;

    Vector<Double> Xprevious_;
    Vector<Double> Yprevious_;

    Matrix<Double> preisachWeights_;
  };


} //end of namespace


#endif

