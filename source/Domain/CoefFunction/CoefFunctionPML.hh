// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     CoefFunctionPML.hh
 *       \brief    This file contains the definitions of the class DampFunction which is used by
 *                 the PML classes, the base PML class CoefFunctionPMLBase, the classic/Cartesian 
 *                 CoefFunctionPML and the (also Cartesian) CoefFunctionShiftedPML
 *
 *       \date     Mar 12, 2012, edit by pheidegger in Mar 2023
 *       \author   ahueppe, pheidegger
 */
//================================================================================================

#ifndef COEFFUNCTIONPML_HH_
#define COEFFUNCTIONPML_HH_

#define _USE_MATH_DEFINES
#include <cmath>

#include "CoefFunction.hh"
#include "Domain/Mesh/Grid.hh"
#include "Utils/mathParser/mathParser.hh"
#include "boost/bind/bind.hpp"
#include "boost/lexical_cast.hpp"
#include "Domain/Mesh/Grid.hh"
#include "Domain/CoordinateSystems/CoordSystem.hh"

namespace CoupledField{

class DampFunction{
  //================================================================================================
  // DampFunction
  //================================================================================================
  //! Class to define damping functions used within the coefFunction*PML* classes
  public:
    typedef enum{ NO_TYPE, CONSTANT, INVERSE_DIST, QUADRATIC, SMOOTH, TANGENS, RATIONAL, EXPONENTIAL, POLY_DIRECT, POLY_INVERSE } DampingType;
    static Enum<DampingType> DampingTypeEnum;

    DampFunction(){
      ReflectionCoefficient = 1e-3;
      functionType = NO_TYPE;
      DampFactor = 1.0;
      constFactor =1.0;
    }

    virtual ~DampFunction() {};

    //! returns the damping function evaluated at position 'pos' in a PML layer with 
    //! thickness 'thickness'
    virtual Double ComputeFactor(Double pos, Double thickness)=0;

    //! returns the integral of the damping function from 0 to the position 'pos'
    virtual Double ComputeIntegralFactor(Double pos, Double thickness) {
      EXCEPTION("ComputeIntegralFactor not overridden by child class " << functionType << "!");
    };

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
      functionType = CONSTANT;
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

    Double ComputeIntegralFactor(Double pos, Double thickness) {
      Double val = 0.0;
      Double a = -log(thickness);
      Double b;
      if(thickness == pos) {
        b = 1e10;//-log(1e-10); // deal with singularity
      } else {
        b = -log(thickness-pos);
      }
      val = constFactor * (b-a);
      return val*DampFactor; // why would I need this DampFactor that is a constant 1???
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
      value *= ( (pos / thickness) - ((sin(2*M_PI*pos / thickness)/(8*atan(1.0))) ) );
      return value*DampFactor;
    }

  };

  class DampFunctionTangens : public DampFunction{
  public:
    DampFunctionTangens( ) : DampFunction(){
      constFactor = 2.0/M_PI;
      functionType = TANGENS;
    }
    // The factor is the derivative of the mapping function
    Double ComputeFactor(Double z, Double sos){
      //Double z = pos/thickness;
      //Double x = DampFactor*tan(z*M_PI/2.0);
      //return 2.0*DampFactor/(M_PI*(DampFactor*DampFactor+x*x));
      Double L = DampFactor*sos; // L corresponds to kappa in the IML paper
      Double c = cos(z/constFactor);
      return c*c/L; // same but possibly faster
    }

  };

  class DampFunctionRational : public DampFunction{
  public:
    DampFunctionRational( ) : DampFunction(){
      constFactor = 1.0;
      functionType = RATIONAL;
    }

    Double ComputeFactor(Double z, Double sos){
      //Double z = pos/thickness; // coordinate in layer
      Double L = DampFactor*sos; // L corresponds to kappa in the IML paper
      Double x = z*L/(1.0 - z); // x-coordinate
      return L/((x+L)*(x+L)); // dz/dx=eta(x)=eta(x(z))
    }

  };

  class DampFunctionExponential : public DampFunction{
  public:
    DampFunctionExponential( ) : DampFunction(){
      constFactor = 1.0;
      functionType = EXPONENTIAL;
    }

    Double ComputeFactor(Double z, Double sos){
      //Double z = pos/thickness; // local coordinate in layer [0,1]
      Double L = DampFactor*sos;
      return (1-z)/L;
    }

  };

  class DampFunctionPolyDirect : public DampFunction
  {
  private:
    UInt power_;

  public:
    DampFunctionPolyDirect(UInt power) : DampFunction()
    {
      power_ = power;
      constFactor = -0.5*(power_ + 1)*log(ReflectionCoefficient);
      functionType = POLY_DIRECT;
    }

    Double ComputeFactor(Double pos, Double thickness)
    {
      Double value = pow(pos/thickness, power_);
      return DampFactor*value;
    }
  };

  class DampFunctionPolyInverse : public DampFunction
  {
  private:
    UInt power_;

  public:
    DampFunctionPolyInverse(UInt power) : DampFunction()
    {
      power_ = power;
      constFactor = -0.5*(power_ + 1)*log(ReflectionCoefficient);
      functionType = POLY_INVERSE;
    }

    Double ComputeFactor(Double pos, Double thickness)
    {
      Double value = pow(1.0 - pos/thickness, power_);
      return DampFactor*value;
    }
  };

  //================================================================================================
  // PML Base Class
  //================================================================================================
  /* Base class to define a CoefFunction that triggers the computation of the 
  /* damping functions and creates the tensors/vectors/scalar values for the 
  /* PML damping
  */
  //! Enumeration data type describing formulations of PML
  enum PMLFormulType { CLASSIC, SHIFTED, CURVILINEAR };

  template<typename T>
  class CoefFunctionPMLBase : public CoefFunction{
  public:
    //! base constructor
    CoefFunctionPMLBase(PtrParamNode pmlDef, PtrCoefFct matCoef,
                    shared_ptr<EntityList> EntList,
                    StdVector<RegionIdType> pdeDomains,
                    CoefFunction::CoefDimType dimType);

    //! destructor
    virtual ~CoefFunctionPMLBase();

    //! get object name as a string
    string GetName() const { return name_; }

    //! return name in string format...
    std::string ToString() const {
        std::string out = this->name_;
        return out;
    };

    //! Return complex-valued tensor at integration point.
    virtual void GetTensor(Matrix<Complex>& tensor, const LocPointMapped& lpm ) {
      EXCEPTION("CoefFunctionPMLBase::GetTensor() not overwritten by child class.");
    }; 

    //! Return a real-valued tensor at integration point.
    virtual void GetTensor(Matrix<Double>& tensor, const LocPointMapped& lpm ) {
      EXCEPTION("CoefFunctionPMLBase::GetTensor() not overwritten by child class.");
    };

    //! Return complex-valued vector at integration point.
    virtual void GetVector(Vector<Complex>& vec, const LocPointMapped& lpm )  {
      EXCEPTION("CoefFunctionPMLBase::GetVector() not overwritten by child class.");
    };

    //! Return real-valued vector at integration point.
    virtual void GetVector(Vector<Double>& vec, const LocPointMapped& lpm )  {
      EXCEPTION("CoefFunctionPMLBase::GetVector() not overwritten by child class.");
    };

    //! Return complex-valued scalar at integration point.
    virtual void GetScalar(Complex& val, const LocPointMapped& lpm ) {
      EXCEPTION("CoefFunctionPMLBase::GetScalar() not overwritten by child class.");
    };

    //! Return real-valued scalar at integration point. This is not implemented as the 
    //! Jakobi determinant will always be complex-valued here
    virtual void GetScalar(Double& val, const LocPointMapped& lpm )  {
      EXCEPTION("CoefFunctionPMLBase::GetScalar() not overwritten by child class.");
    };

    //! \copydoc CoefFunction::GetVecSize
    virtual UInt GetVecSize() const {
      assert(this->dimType_ == CoefFunction::VECTOR );
      return this->dim_;
    }

    //! \copydoc CoefFunction::GetTensorSize
    virtual void GetTensorSize( UInt& numRows, UInt& numCols ) const {
      assert(this->dimType_ == CoefFunction::TENSOR );
      numRows = this->dim_;
      numCols = this->dim_;
    }

    //! disable manually adding an entity list, as the entity list is set in the constructor
    void AddEntityList(shared_ptr<EntityList>){
      EXCEPTION("Add Entities may not be called in PML CoefFunction. Specify the region in the constructor!");
    };

    //! return if this instance is called as complex
    bool IsComplex(){
      return std::is_same<T,Complex>::value;
    };

    //! check for current angular frequency
    void UpdateOmega();

    //! dummy function to read the PML data. Must be implemented in the child classes 
    virtual void ReadDataPML(PtrParamNode pmlDef,StdVector<RegionIdType> pdeDomains) {
      EXCEPTION( "CoefFunction::ReadDataPML not overwritten by " << GetName());
    };

    //! Set the type of damping function in the damping-function object 
    void CreateDampFunction();
  protected:
      //! name of the CoefFunctionPML
      std::string name_;
      //! PML formulation type
      PMLFormulType formulationType_;
      //! Support of the CoefFunction. Only needed for grid/solution results
      StdVector<shared_ptr<EntityList>> entities_;
      //! pointer to an instance of the DampFunction class
      shared_ptr<DampFunction> dampFunction_;
      //! type of the damping function
      DampFunction::DampingType pmlType_;
      //! speed of wave
      PtrCoefFct matCoef_;
      //! Pointer to math parser instance
      MathParser* mp_;
      //! Handle for expression
      unsigned int mHandle_;
      //! storing the current frequency
      Double omega_;
      //! dimension of the problem
      UInt dim_;
      //! order of the coeff function: 0->scalar, 1->vector, 2->tensor
      UInt orderCoefFct_;
  };

  //================================================================================================
  // Classic (=Cartesian) PML
  //================================================================================================
  /* This class represents the original 'classic' PML coefFunction, 
  /* which computes the damping vectors/scalars for the Cartesian PML
  */
  template<typename T>
  class CoefFunctionPML : public CoefFunctionPMLBase<T> {

  public:
    CoefFunctionPML(PtrParamNode pmlDef, PtrCoefFct matCoef,
                    shared_ptr<EntityList> EntList,
                    StdVector<RegionIdType> pdeDomains,
                    CoefFunction::CoefDimType dimType);

    virtual ~CoefFunctionPML();

    //! Return real-valued tensor at integration point
    virtual void GetTensor(Matrix<Complex>& tensor,
                  const LocPointMapped& lpm ) override; 

    //! Return real-valued tensor at integration point
    virtual void GetTensor(Matrix<Double>& tensor,
                  const LocPointMapped& lpm ) override;

    //! Return real-valued vector at integration point
    virtual void GetVector(Vector<Complex>& vec,
                  const LocPointMapped& lpm ) override;

    //! Return real-valued vector at integration point
    virtual void GetVector(Vector<Double>& vec,
                  const LocPointMapped& lpm ) override;

    //! Return real-valued scalar at integration point
    // this is little bit of a hack,
    // seeing that the jacobian is transformed according to the changed
    // derivatives, we pass this function as a scalar function to the bilinearform
    // an transform the jacobian with it....
    virtual void GetScalar(Double& val,
                  const LocPointMapped& lpm ) override;

    //! Return cpmplex-valued scalar at integration point
    virtual void GetScalar(Complex& val,
                  const LocPointMapped& lpm ) override;

    //! use the base functions
    using CoefFunctionPMLBase<T>::UpdateOmega;
    using CoefFunctionPMLBase<T>::CreateDampFunction;

  protected:
    void ReadDataPML(PtrParamNode pmlDef,StdVector<RegionIdType> pdeDomains) override;

    void GuessLayerData(StdVector<RegionIdType> pdeDomains);

    void GetThicknessAtPoint(Double& thickness,Double& position, LocPointMapped lpm,UInt dir);

    Matrix<Double> innerMinMaxComp_;
    Matrix<Double> outerMinMaxComp_;

  };


  //================================================================================================
  // Shifted PML
  //================================================================================================
  template<typename T>
  class CoefFunctionShiftedPML : public CoefFunctionPML<T>
  {

  public:

  CoefFunctionShiftedPML(PtrParamNode pmlDef, PtrCoefFct matCoef, shared_ptr<EntityList> EntList,
                         StdVector<RegionIdType> pdeDomains, CoefFunction::CoefDimType dimType);

    virtual ~CoefFunctionShiftedPML();

    //! \copydoc CoeffFunctionPML::GetTensor
    virtual void GetTensor(Matrix<Complex>& tensor, const LocPointMapped& lpm) override;

    //! \copydoc CoeffFunctionPML::GetVector
    virtual void GetVector(Vector<Complex>& vector, const LocPointMapped& lpm) override;

    //! \copydoc CoeffFunctionPML::GetScalar
    virtual void GetScalar(Complex& scalar, const LocPointMapped& lpm) override;

    //! use the functions defined in the CoefFunctionPML, e.g. for GetTensor(Matrix<Double>....)
    using CoefFunctionPML<T>::GetTensor;
    using CoefFunctionPML<T>::GetVector;
    using CoefFunctionPML<T>::GetScalar;

  private:
    PtrCoefFct scalingCoef_, shiftCoef_;
    shared_ptr<DampFunction> scalingFunc_, shiftFunc_;
  };
}
#endif /* COEFFUNCTIONPML_HH_ */
