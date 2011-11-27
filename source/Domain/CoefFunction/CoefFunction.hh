// =====================================================================================
// 
//       Filename:  coefFunction.hh
// 
//    Description:  This is the base class for every coefficient
//                  It is used in the Integrator classes to obtain the material tensor
//                  (D-Matrix) as well as in the B-Operator classes to realize e.g. a
//                  change in the derivatives. 
//
//                  Example1: consider the PML
//                  Coefficient function provides a material matrix multiplied by the
//                  complex jacobian. Moreover the same coeffient function is passed to
//                  the operator to provide the daming parameters for the "streched"
//                  derivatives.
//
//                  Example2: Consider geometric non-linear mechanics
//                  The CoefFunction can provide the solution at the current integration
//                  point to enable the b-operator to compute its matrix of derivatives
//
//                  This class is just the interface class right now, the following
//                  structure should be realized:
//                  - coefFunctionAnalytic
//                    This class descibes the quantity by an analytic expression this
//                    includes (non-linear) materials, analytic flow fields, etc.
//                  - coefFunctionGrid
//                    This class provides the Information based on information read in
//                    from a grid. E.g. in the case of (interpolated) aeroacoustic source terms 
//                  - coefFucnctionSol
//                    This function provides the information based on the current
//                    solution of the problem
//
//
// 
//        Version:  1.0
//        Created:  10/19/2011 10:16:34 AM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Andreas Hueppe (AHU), andreas.hueppe@uni-klu.ac.at
//        Company:  Universitaet Klagenfurt
// 
// =====================================================================================

#ifndef COEFFUNCTION_HH
#define COEFFUNCTION_HH

#include "General/Environment.hh"
#include "MatVec/Matrix.hh"
#include "MatVec/Vector.hh"
#include "Domain/Domain.hh"
#include "Domain/ElemMapping/Elem.hh"
#include "Domain/ElemMapping/ElemShapeMap.hh"
#include "Domain/ElemMapping/EntityLists.hh"
namespace CoupledField{


// forward class declarations
class CoordSystem;

class CoefFunction{
public:

  // ========================
  //  ENUMS / TYPEDEFS
  // ========================
  //@{ \name Public type definitions

  //! Dimension of coefficient functions
  typedef enum{ 
    NO_DIM,   /*!< Uninitialized */
    SCALAR,   /*!< 0-dimensional scalar entry */
    VECTOR,   /*!< 1-dimensional vector entry */
    TENSOR    /*!< 2-dimensional tensor entry */ 
  } CoefDimType;
  
  //! Dependency of coefficient function
  typedef enum{ 
    CONST,         /*!< No dependency on space or time */
    TIMEFREQ,      /*!< Only depending on time / frequency, not space */
    SPACE,         /*!< Only spatial dependency */
    SPACE_TIMEFREQ /*!< Spatial and time/frequency dependency */
  } CoefDependType;
  
  // ========================
  //  FACTORY METHODS
  // ========================
  //@{ \name Factory Methods

  //! Generate scalar-valued coefficient function
  static shared_ptr<CoefFunction> 
  Generate( Global::ComplexPart format, 
            const std::string& realVal, 
            const std::string& imagVal = "" );

  //! Generate vector-valued coefficient function
  static shared_ptr<CoefFunction> 
  Generate( Global::ComplexPart format, 
            const StdVector<std::string>& realVal, 
            const StdVector<std::string>& imagVal = StdVector<std::string>() );


  //! Generate tensor-valued coefficient function
  static shared_ptr<CoefFunction> 
  Generate( Global::ComplexPart format,
            UInt numRows, UInt numCols,
            const StdVector<std::string>& realVal,
            const StdVector<std::string>& imagVal = StdVector<std::string>() );
  //@}

  // ========================
  //  CONSTRUCTOR / DESTRUCTOR
  // ========================
  //@{ \name Constructor / Destructor

  CoefFunction(){
    dimType_ = NO_DIM;
    dependType_ = CONST;
    //set the default coordinate system
    // can be reset later by user...
    this->coordSys_ = domain->GetCoordSystem();
  }

  ~CoefFunction(){
    ;
  }
  //@}
  
  // ========================
  //  ACCESS METHODS
  // ========================
  //@{ \name Access Methods

  virtual void GetTensor(Matrix<Double>& CoefMat, 
                         const LocPointMapped& lpm ) {
    EXCEPTION("CoefFunction::GetTensor<Double> called: This may not happen");
  }

  virtual void GetVector(Vector<Double>& CoefMat, 
                         const LocPointMapped& lpm ) {
    EXCEPTION("CoefFunction::GetVector<Double> called: This may not happen");
  }

  virtual void GetScalar(Double& CoefMat, 
                         const LocPointMapped& lpm ) {
    EXCEPTION("CoefFunction::GetScalar<Double> called: This may not happen");
  }

  virtual void GetTensor(Matrix<Complex>& CoefMat, 
                         const LocPointMapped& lpm ) {
    EXCEPTION("CoefFunction::GetTensor<Complex> called: This may not happen");
  }

  virtual void GetVector(Vector<Complex>& CoefMat, 
                         const LocPointMapped& lpm ) {
    EXCEPTION("CoefFunction::GetVector<Complex> called: This may not happen");
  }

  virtual void GetScalar(Complex& CoefMat, 
                         const LocPointMapped& lpm ) {
    EXCEPTION("CoefFunction::GetScalar<Complex> called: This may not happen");
  }
  //@}


  void SetCoordinateSystem(CoordSystem* cSys){
    coordSys_ = cSys;
  }

  virtual CoefDimType GetDimType(){
    return dimType_;
  }

  virtual void AddEntities(shared_ptr<EntityList> ent){
    if(!entities_.Contains(ent)){
      entities_.Push_back(ent);
    }else{
      WARN("entity list " << ent->GetName() << " already contained in CoefFunction")
    }
  }

  virtual void UpdateCoefs(Double tfParam){
    WARN(__FILE__ << __LINE__ << " Update Coefs has no functionality yet");
  }

  //! Dump coefficient function to string 
  virtual std::string ToString() const {
    EXCEPTION("CoefFuncion: ToString() not properly overwritten");
  }


protected:

  // ========================
  //  HELPER METHODS
  // ========================
  //@{ \name Helper methods
  
  //! Returns true, if expression depends on time / freq
  static bool ExprDependsOnTimeFreq(const std::string& expr);
  
  //! Returns true, if expression depends on space
  static bool ExprDependsOnSpace(const std::string& expr);
  //@}
  
  //TODO: CHANGE THIS TO SHARED POINTER
  // i.e. change the daomin to hold a shared_ptr to the coordinate systems!
  CoordSystem* coordSys_;

  //! Dimension of coefficient function (scalar, vector, tensor)
  CoefDimType dimType_;
  
  //! Dependency type of the coefficient function
  CoefDependType dependType_;

  //! Support of the CoefFunction. Only needed for grid/solution results
  StdVector<shared_ptr<EntityList> > entities_;
};

}
#endif
