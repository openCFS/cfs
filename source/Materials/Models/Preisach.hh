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

    //!
    virtual ~Preisach();

    //! destructor
    Double computeValue(Double& xVal, Integer idxElem, bool overwrite = true);

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
    Double normalizeInput(Double xInput);

    //! normalizes the output to Ysaturated_
    Double normalizeOutput(Double xInput);

    std::string runtimeToString(){
      return "No runtime information available for Scalar Preisach model";
    };

    void setFlags(UInt performanceFlag, UInt textOutputFlag, UInt mappingFlag){
      ;
    };

  protected:

    //! computes  the everett function (area-integration for x1, x2)
    Double everettPixel(Double x1, Double x2);

  private:

    Double Xsaturated_; //! saturation value for  input
    Double YSaturated_; //! saturation value for output

    bool isVirgin_; //! yes, if starting at zero

    Vector<Double> preisachSum_; //! output value of Preisach operator

    Vector<Double>* strings_; //! irreduceable minima and maxima
    Vector<Double>* helpStrings_; //! help array for string_

    Vector<UInt> StringLength_; //! number of irreduceable minima and maxima
    UInt maxStringLength_; //! maximum allowd length for 

    Matrix<Double> preisachWeights_; //! presach weight function

    Double eps_; //! accuracy parameter
  };


} //end of namespace


#endif

