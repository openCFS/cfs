// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     CoefFunctionPML.hh
 *       \brief    <Description>
 *
 *       \date     Mar 12, 2012
 *       \author   ahueppe
 */
//================================================================================================

#ifndef COEFFUNCTIONPML_HH_
#define COEFFUNCTIONPML_HH_

#include "CoefFunction.hh"
#include <boost/tr1/type_traits.hpp>

namespace CoupledField{

class DampFunction{

public:
  typedef enum{ NO_TYPE, CONST, INVERSE_DIST, QUADRATIC, SMOOTH } DampingType;
  static Enum<DampingType> DampingTypeEnum;

  DampFunction(){
    ReflectionCoefficient = 1e-3;
    functionType = NO_TYPE;
    DampFactor = 1.0;
  }

  virtual ~DampFunction() {};

  virtual Double ComputeFactor(Double pos, Double thickness)=0;


  DampingType GetType(){return functionType;}

  Double DampFactor;

protected:
  Double constFactor;
  Double ReflectionCoefficient;
  DampingType functionType;


};

class DampFunctionConst : public DampFunction{
public:
  DampFunctionConst(Double SpeedOfSound) : DampFunction(){
    constFactor = -1.0 * SpeedOfSound * log(ReflectionCoefficient)/2.0;
    functionType = CONST;
  }

  Double ComputeFactor(Double pos, Double thickness){
    return DampFactor * constFactor / thickness;
  }

};

class DampFunctionQuad : public DampFunction{
public:
  DampFunctionQuad(Double SpeedOfSound) : DampFunction(){
    constFactor = -3.0 * SpeedOfSound * log(ReflectionCoefficient)/2.0;
    functionType = QUADRATIC;
  }

  Double ComputeFactor(Double pos, Double thickness){
    return DampFactor*constFactor * pos * pos / (thickness * thickness*thickness);
  }

};

class DampFunctionInvDist : public DampFunction{
public:
  DampFunctionInvDist(Double SpeedOfSound) : DampFunction(){
    constFactor = SpeedOfSound;
    functionType = INVERSE_DIST;
  }

  Double ComputeFactor(Double pos, Double thickness){
    //check for sigular position
    Double val = 0.0;
    if(thickness == pos){
      //ok we just take a high value....
      val = constFactor * (1.0 / 1e-10);
    }else{
      val = constFactor * (1.0 / (thickness - pos));
    }
    return val*DampFactor;
  }

};

class DampFunctionSmooth : public DampFunction{
public:
  DampFunctionSmooth(Double SpeedOfSound) : DampFunction(){
    constFactor = SpeedOfSound * log(1.0/ReflectionCoefficient);
    functionType = SMOOTH;
  }

  Double ComputeFactor(Double pos, Double thickness){
    Double value = constFactor/thickness;
    value *= ( (pos / thickness) - ((sin(2*M_PI*pos / thickness)/(2*M_PI)) ) );
    return value*DampFactor;
  }

};



template<typename T>
class CoefFunctionPML : public CoefFunction{

public:


  CoefFunctionPML(PtrParamNode pmlDef, PtrCoefFct speedOfSound,
                  shared_ptr<EntityList> EntList,
                  StdVector<RegionIdType> pdeDomains,
                  bool isVector );

  ~CoefFunctionPML();

  //! Return real-valued tensor at integration point
  void GetTensor(Matrix<Complex>& tensor,
                 const LocPointMapped& lpm );

  //! Return real-valued tensor at integration point
  void GetTensor(Matrix<Double>& tensor,
                 const LocPointMapped& lpm );

  //! Return real-valued vector at integration point
  void GetVector(Vector<Complex>& vec,
                 const LocPointMapped& lpm ) ;

  //! Return real-valued vector at integration point
  void GetVector(Vector<Double>& vec,
                 const LocPointMapped& lpm ) ;


 //! Return real-valued scalar at integration point
  // this is little bit of a hack,
  // seeing that the jacobian is transformed according to the changed
  // derivatives, we pass this function as a scalar function to the bilinearform
  // an transform the jacobian with it....
 void GetScalar(Double& val,
                const LocPointMapped& lpm ) ;

 //! Return cpmplex-valued scalar at integration point
 void GetScalar(Complex& val,
                const LocPointMapped& lpm ) ;

 //! \copydoc CoefFunction::GetVecSize
 UInt GetVecSize() const {
   assert(this->dimType_ == VECTOR );
   return dim_;
 }
 

  void AddEntities(shared_ptr<EntityList>){
    EXCEPTION("Add Entities may not be called in PML CoefFunction. Specify the region in the constructor!");
  }

  virtual bool IsComplex(){
    return std::tr1::is_same<T,Complex>::value;
  }

protected:
    //void ComputeDampingFactor( Vector<Double>& factors,
    //                           const LocPointMapped& lpm,
    //                           UInt dir);
    //
    //void ComputeDampingFactor(Vector<Complex>& factors,
    //                          const LocPointMapped& lpm,
    //                          UInt dir);
    //
    //! Call-back method for re-calculation
    void UpdateOmega(){
      omega_ = this->mp_->Eval(mHandle_) * 2 * M_PI;
    }

private:
    void SetPosPML(Matrix<Double> & inner,
                   Matrix<Double> & outer);

    void ReadDataPML(PtrParamNode pmlDef,StdVector<RegionIdType> pdeDomains);

    void CreateDampFunction();

    void GuessLayerData(StdVector<RegionIdType> pdeDomains);

    void GetThicknessAtPoint(Double& thickness,Double& position, LocPointMapped lpm,UInt dir);


    Matrix<Double> innerMinMaxComp_;
    Matrix<Double> outerMinMaxComp_;

    shared_ptr<DampFunction> dampFunction_;

    DampFunction::DampingType pmlType_;

    //! should be coefFunction in future...
    PtrCoefFct speedOfSound_;

    //! Pointer to math parser instance
    MathParser* mp_;

    //! Handle for expression
    MathParser::HandleType mHandle_;

    //! storing the current frequency
    Double omega_;
    
    //! dimension of the problem
    UInt dim_;

    //! flag, if PML coefficient functions describes the vector 
    bool isVector_;

};

}
#endif /* COEFFUNCTIONPML_HH_ */
