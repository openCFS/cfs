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

#define _USE_MATH_DEFINES
#include <cmath>

#include "CoefFunction.hh"
#include "Domain/Mesh/Grid.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include "FeBasis/H1/H1Elems.hh"
//#include "FeBasis/FeSpace.hh"

namespace CoupledField{

  // forward class declaration, needed in CoefFunctionCurvilinearPML
  //class GridCFS;

  //================================================================================================
  // DampFunction
  //================================================================================================
  /* Class to define damping functions used within the coefFunction*PML* classes
  */
  class DampFunction{

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
    virtual Double ComputeFactor(Double pos, Double thickness) {
      EXCEPTION("ComputeFactor not overridden by child class " << functionType << "!");
    };

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

    Double ComputeFactor(Double pos, Double thickness) {
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
    CoefFunctionPMLBase(PtrParamNode pmlDef, PtrCoefFct speedOfSound,
                    shared_ptr<EntityList> EntList,
                    StdVector<RegionIdType> pdeDomains);

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
      //! speed of sound
      PtrCoefFct speedOfSound_;
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
    CoefFunctionPML(PtrParamNode pmlDef, PtrCoefFct speedOfSound,
                    shared_ptr<EntityList> EntList,
                    StdVector<RegionIdType> pdeDomains,
                    bool isVector );

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
    void ReadDataPML(PtrParamNode pmlDef,StdVector<RegionIdType> pdeDomains);

    void GuessLayerData(StdVector<RegionIdType> pdeDomains);

    void GetThicknessAtPoint(Double& thickness,Double& position, LocPointMapped lpm,UInt dir);

    Matrix<Double> innerMinMaxComp_;
    Matrix<Double> outerMinMaxComp_;

    //! flag, if PML coefficient functions describes the vector 
    bool isVector_;
  };


  //================================================================================================
  // Shifted PML
  //================================================================================================
  template<typename T>
  class CoefFunctionShiftedPML : public CoefFunctionPML<T>
  {

  public:

    CoefFunctionShiftedPML(PtrParamNode pmlDef, PtrCoefFct speedOfSound, shared_ptr<EntityList> EntList,
                          StdVector<RegionIdType> pdeDomains, bool isVector);

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


  //================================================================================================
  // Curvilinear PML
  //================================================================================================
  /* A curvilinear PML that is supposed to function in arbitrarily-shaped, convex 
  * domains. 
  * The CoefFunctionCurvilinearPML currently works in combination with an automatic layer generation
  * and the triggers the computation of geometry of the layer (normal vectors, principal-direction 
  * vectors and principal curvatures) in the grid class. 
  * The computed nodal data is interpolated to the integration points and the scalars/tensors
  * are assembled directly there.
  */

  //! the type of output the coefFunction is computing (currently only for curvilinearPML):
  //! TENSOR (tensor output) ... the whole factorized PML matrix J^-1 = A^T (I-D)^-1 A
  //! DETERMINANT (scalar output) ... the determinant of the PML matrix
  //! DAMP_FACTOR (vector output) ... the eigenvalues of the PML matrix (entries of the 'inner' diagonal matrix outputted as vector)
  //! NORM_VEC (vector output) ... normal vecors within the PML layer, computed on nodes
  //! MIN_PRINC_VEC (vector output) ... minimum principal direction vector within the layer, computed on nodes
  //! MIN_PRINC_CURV (scalar output) ... minimum principal curvature within the layer, computed on nodes
  //! MAX_PRINC_VEC (vector output) ... maximum principal direction vector within the layer, computed on nodes (unused in 2D)
  //! MAX_PRINC_CURV (scalar output) ... maximum principal curvature within the layer, computed on nodes (unused in 2D)
  //! DISTANCE (scalar output) ... the distance of a point to its closest point on the PML interface
  enum OutputType { TENSOR, DETERMINANT, DAMP_FACTOR, NORM_VEC, MIN_PRINC_VEC, MIN_PRINC_CURV, MAX_PRINC_VEC, MAX_PRINC_CURV, DISTANCE};
  template<typename T>
  class CoefFunctionCurvilinearPML : public CoefFunctionPMLBase<T> {
  public:
    //! constructor
    //! \param pmlDef (in) param node that holds the PML parameters defined in the XML
    //! \param speedOfSound (in) coef function with the speed of sound (gets assigned to member speedOfSound in the base class)
    //! \param EntList (in) entity list containing all entities of the assigned PML region
    //! \param pdeDomains (in) StdVector holding all domains of the current PDE
    //! \param outputType (in) enum to set which output should be computed
    CoefFunctionCurvilinearPML(PtrParamNode pmlDef, PtrCoefFct speedOfSound, shared_ptr<EntityList> EntList,
                          StdVector<RegionIdType> pdeDomains, OutputType outputType);

    //! destructor
    virtual ~CoefFunctionCurvilinearPML();

    //! Return the complex-valued tensor of the PML at integration point. The Tensor will be used in 
    //! TensorScaledGradientOperator() to compute the scaled derivative (for AcousticPDE).
    //! This function gets only called when the outputType_ is TENSOR
    //! \param tensor (out) the tensor holding the complex-valued inverse Jakobi matrix of
    //!                     the curvilinear coordinate stretch
    //! \param lpm (in) the local point for which the tensor is returned
    virtual void GetTensor(Matrix<Complex>& tensor, const LocPointMapped& lpm ) override; 

    //! Return a real-valued tensor at integration point. This is not implemented as the 
    //! coordinate stretching operation is generally complex-valued.
    virtual void GetTensor(Matrix<Double>& tensor, const LocPointMapped& lpm ) override {
      EXCEPTION("CoefFunctionCurvilinearPML::GetTensor() is only implemented for Complex values.");
    }

    //! Return complex-valued scalar at integration point. Returns one of the four vector-valued outputTypes, 
    //! depending on which type is set.
    //! This function gets only called when the outputType_ is DAMP_FACTOR, NORM_VEC, MIN_PRINC_VEC, or MAX_PRINC_VEC.
    virtual void GetVector(Vector<Complex>& vec, const LocPointMapped& lpm ) override;

    //! Return real-valued scalar at integration point. Returns one of the three vector-valued outputTypes, 
    //! depending on which type is set.
    //! This function gets only called when the outputType_ is NORM_VEC, MIN_PRINC_VEC, or MAX_PRINC_VEC.
    virtual void GetVector(Vector<Double>& vec, const LocPointMapped& lpm ) override;

    //! Return complex-valued scalar at integration point. Returns either the determinant of the PML Jacobi matrix,
    //! or one of the other three scalar-valued outputTypes, depending on which type is set.
    //! This function gets only called when the outputType_ is DETERMINANT, DISTANCE, MIN_PRINC_CURV, or MAX_PRINC_CURV
    //! \param val (out) the scalar holding the complex-valued Jakobi determinant of the 
    //!                  curvilinear coordinate stretch
    //! \param lpm (in) the local point for which the scalar is returned
    virtual void GetScalar(Complex& val, const LocPointMapped& lpm ) override;

    //! Return real-valued scalar at integration point. Returns one of the three scalar-valued outputTypes, 
    //! depending on which type is set.
    //! This function gets only called when the outputType_ is DISTANCE, MIN_PRINC_CURV, or MAX_PRINC_CURV
    virtual void GetScalar(Double& val, const LocPointMapped& lpm ) override;

    //! use the base functions
    using CoefFunctionPMLBase<T>::UpdateOmega;
    using CoefFunctionPMLBase<T>::CreateDampFunction;

  private:
    //! read data from the paramNode and store
    void ReadDataPML(PtrParamNode pmlDef);

    //! check if autoLayerGeneration is specified and if it fits to the current PML region.
    //! If true, sets the layerGenNode_, numLayers_, elemHeight_, numSurfNodes_, and layerThickness_.
    void CheckForLayerGenerationNode(PtrParamNode pmlDef);

    //! check if invalid parameters are set and warn the user
    void CheckForInvalidParams(PtrParamNode pmlDef);

    //! compute distances of every node in the layer to its closest point on the PML interface.
    //! Stores the values in thicknessOnNodes_.
    void GetThicknessOnNodes();

    //! returns the position of a given nodeId in the NodeGeometry struct (which is supposed to 
    //! equal its index in the PML region). Exploits the knowledge of the auto layer generation params
    //! so speed up the search.
    UInt GetIdxByNodeId(const UInt& nodeId) const;

    //! Interpolates the nodal quantities to the given local point mapped. 
    //! sets: n_, tmin_, tmax_, kmin_, kmax_, dist_, sos_, dampFunc_, and intDampFunc_; (for Tensor in 3D)
    //! If the CoefFunction is used as Scalar, only kmin_, kmax_, sos_, dist_, dampFunc_, and intDampFunc are set.
    //! In 2D, tmin_ and kmax_ are used for tangential vector and curvature, tmax_ and kmax_ are not set.

    //! Interpolates the nodal quantities to the given local point mapped. 
    //! Sets: n_, tmin_, tmax_, kmin_, kmax_, dist_, sos_, dampFunc_, and intDampFunc_;
    void GetTensorParams(const LocPointMapped& lpm);

    //! Interpolates the nodal quantities to the given local point mapped. 
    //! Sets: kmin_, kmax_, dist_, sos_, dampFunc_, and intDampFunc_;
    void GetDeterminantParams(const LocPointMapped& lpm);

    //! Interpolates the nodal quantities to the given local point mapped. 
    //! Sets: kmin_, kmax_, dist_, sos_, dampFunc_, and intDampFunc_;
    void GetDampFactorParams(const LocPointMapped& lpm);

    //! Interpolates the nodal quantities to the given local point mapped. 
    //! Sets: n_;
    void GetNormVecParams(const LocPointMapped& lpm);

    //! Interpolates the nodal quantities to the given local point mapped. 
    //! Sets: tmin_;
    void GetMinPrincVecParams(const LocPointMapped& lpm);

    //! Interpolates the nodal quantities to the given local point mapped. 
    //! Sets: kmin_;
    void GetMinPrincCurvParams(const LocPointMapped& lpm);

    //! Interpolates the nodal quantities to the given local point mapped. 
    //! Sets: tmax_;
    void GetMaxPrincVecParams(const LocPointMapped& lpm);

    //! Interpolates the nodal quantities to the given local point mapped. 
    //! Sets: kmax_;
    void GetMaxPrincCurvParams(const LocPointMapped& lpm);

    //! Interpolates the nodal quantities to the given local point mapped. 
    //! Sets: dist_;
    void GetDistParams(const LocPointMapped& lpm);

    //! Create identity operators for mapping the nodal values (e.g. geometry data) to local points.
    //! Defines one operator for vectors: vectorMappingOperator_; and one for scalars: scalarMappingOperator_;
    void CreateMappingOperators();

    //! Set the dimType_ of the coefFunction according to the assigned outputType_
    void SetDimType();

    void ComputeTensorOnNode(const UInt& nodeIdx);

    //StdVector<Matrix<Complex>> tensorsOnNodes_;

    //! member to set which quantity should be outputted/computed by the coefFunction
    OutputType outputType_;

    //! pointer to the current grid
    Grid* grid_;

    //! pointer to the layer generation node that keeps auto-mesh-generation parameters 
    //! (elemHeight and numLayers). Is empty if not specified.
    PtrParamNode layerGenNode_;

    //! the volume region to act on
    RegionIdType volRegion_;

    //! pointer to the stored geometry of every node in the PML layer
    shared_ptr<Grid::NodeGeometry> ptrNodeGeom_;

    //! vector containing the relative distance of every node to its closest 
    //! point on the PML interface
    StdVector<Double> thicknessOnNodes_;

    //! BOperator to map solutions to arbitrary points. Used is an identity operator,
    //! defined for interpolating Double-valued vectors and scalars. Tensors must be 
    //! stacked into one dimension 
    shared_ptr<BaseBOperator > tensorMappingOperator_;
    shared_ptr<BaseBOperator > vectorMappingOperator_;
    shared_ptr<BaseBOperator > scalarMappingOperator_;
  
    //! Layer generation parameters
    Double numLayers_;      //number of generated surface regions within the layer
    Double elemHeight_;     //height of a generated element
    Double numSurfNodes_;   //number of nodes on every surface
    Double layerThickness_; //total thickness of the PML layer

    //! variables to store quantities on lpm
    Vector<Double> n_;       //normal vector
    Vector<Double> tmin_;    //min principal vector
    Vector<Double> tmax_;    //max principal vector
    Double kmin_;            //min principal curvature
    Double kmax_;            //max principal curvature
    Double dist_;            //distance to interface
    Double sos_;             //speed of sound
    Double dampFunc_;        //damping function
    Double intDampFunc_;     //integral over damping function

    //! store geometry of a current node
    Matrix<Double> nMat_;
    Matrix<Double> tminMat_;
    Matrix<Double> tmaxMat_;
  };
}
#endif /* COEFFUNCTIONPML_HH_ */
