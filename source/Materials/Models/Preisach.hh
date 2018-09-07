#ifndef FILE_PREISACH_2004
#define FILE_PREISACH_2004

#include "Hysteresis.hh"

#include <list>

#include "MatVec/Vector.hh"
#include "MatVec/Matrix.hh"


namespace CoupledField {

  class Preisach : public Hysteresis
  {
  public:
    //! constructor
    Preisach(Integer numElem, Double xSat, Double ysat, 
	     Matrix<Double>& preisachWeight, bool isVirgin);

    //! destructor
    virtual ~Preisach();

    //! actually never used
    //Double computeValue(Double& xVal, Integer idxElem, bool overwrite = true);
    
    Double computeInputAndUpdate(Double Yin, Double eps_mu, Integer idx, bool overwrite = true);
    
    //!computes for xVal a new output value and deletion rules are applied
    Double computeValueAndUpdate(Double xVal, Integer idxElem, 
                                 bool overwrite = true);

    //! returns the current output of the hyst-operator for element idxElem
    Double getValue(  Integer idxElem);

    //! returns the number of stored minima and maxima
    UInt getStringLength( Integer idxElem ) {
      return  StringLength_[ idxElem];
    };

    //! updates the list of minima and maxima due to new input
    Double updateMinMaxList(Double newX, Integer idxElem, 
                          bool overwrite);

    //!
    void SetTimeStepVal(Double dt) 
    {;};

    //! normalizes the input to Xsaturated_
    Double normalizeAndClipInput(Double xInput);

    std::string runtimeToString(){
      return "No runtime information available for Scalar Preisach model";
    };

    void setFlags(UInt performanceFlag){
      ;
    };

  protected:

    //! computes  the everett function (area-integration for x1, x2)
    Double everettPixel(Double x1, Double x2);

    Double bisect(Double dY,Double xMin,Double xMax, Double xFixed, Double eps_mu, Double tol);
    
  private:

    Double XSaturated_; //! saturation value for  input
    Double YSaturated_; //! saturation value for output

    bool isVirgin_; //! yes, if starting at zero

    Vector<Double> preisachSum_; //! output value of Preisach operator

    Vector<Double>* strings_; //! irreduceable minima and maxima
    Vector<Double>* helpStrings_; //! help array for string_
    Vector<Integer>* minmaxtype_; // stores for each entry of the min/max list if it is a minimum or not
    // -1 = min; 1 = max; 0 = initial state
    Vector<Double>* evaluatedEverettPixel_; // stores to each entry of the min/max list
                                            // the corresponding everett pixel (with sign)
    // i.e. everett(-strings[0],strings[0]), everett(strings[0],strings[1]), everett(strings[1],strings[2]), ...
    // note that we need one additional entry for the first min/max

    Vector<UInt> StringLength_; //! number of irreduceable minima and maxima
    UInt maxStringLength_; //! maximum allowd length for 

    Matrix<Double> preisachWeights_; //! preisach weight function

    Double tol_; //! accuracy parameter
    
    // previous input value X and polarization P
    // NOTE:
    //  X and P are normalized; P is clipped to saturation, X must not be clipped
    Vector<Double> previousXval_;
    Vector<Double> previousPval_;
  };


} //end of namespace


#endif

