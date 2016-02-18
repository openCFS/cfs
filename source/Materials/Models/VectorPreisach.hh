#ifndef FILE_VECPREISACH_2016
#define FILE_VECPREISACH_2016

#include "MatVec/Vector.hh"
#include "MatVec/Matrix.hh"
#include "Hysteresis.hh"
//#include "Preisach.hh"


namespace CoupledField {

  class VectorPreisach : public Hysteresis
  {
  public:
    //! constructor
    VectorPreisach(Integer numElem, Double xSat, Double ysat,
	     Matrix<Double>& preisachWeight, Double rotationalResistance , UInt dim, bool isVirgin);

    //!
    virtual ~VectorPreisach();

    //! update rotational Weights
    void updateRotationalWeights(Vector<Double>& Xin, Integer idx);

    //! update switching Operators; has to be called after updateRotationalWeights!
    void updateSwitchingOperators(Vector<Double>& Xin, Integer idx);

    //! destructor
    Vector<Double> computeValue_vec(Vector<Double>& xVal, Integer idxElem);

    //! returns the current output of the hyst-operator for element idxElem
    Vector<Double> getValue_vec(  Integer idxElem);


  private:

    Double Xsaturated_; //! saturation value for  input
    Double YSaturated_; //! saturation value for output

    bool isVirgin_; //! yes, if starting at zero

    Vector<Double>* preisachSum_; //! output value of Preisach operator

    Matrix<Double> preisachWeights_; //! preisach weight function
    Matrix<Double>* rotationalWeightsX_; //! rotational operator
    Matrix<Double>* rotationalWeightsY_; //! rotational operator
    Matrix<Double>* rotationalWeightsZ_; //! rotational operator for 3D
    Matrix<Integer>* switchingStates_; //! current state of the switching operators (evaluation via Everett not possible)

    Double eps_; //! accuracy parameter
    Double rotationalResistance_;

    bool* rotationalWeightsUpdated_;
    bool* switchingStatesUpdated_;

    UInt dim_; //! 2D or 3D

  };


} //end of namespace


#endif

