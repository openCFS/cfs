#ifndef FILE_FLUIDMECHSTIFFINT
#define FILE_FLUIDMECHSTIFFINT

#include <string>

#include "General/defs.hh"
#include "linSurfForm.hh"
#include "fluidMechInt.hh"

namespace CoupledField {
class EntityIterator;
template <class TYPE> class Matrix;
}  // namespace CoupledField

namespace CoupledField
{
/**********************************************************
 * ************* Stabilized NewtonFEM ****************
 **********************************************************
 */

class FluidMechPlaneStiffNewtonInt_UV : public FluidMechInt
{
  public:
    FluidMechPlaneStiffNewtonInt_UV(Double density, \
        Double kinematicViscosity,Matrix<Double> & stabilParams, \
        bool movingMesh, std::string stabilType);
    virtual ~FluidMechPlaneStiffNewtonInt_UV();
    void CalcElementMatrix( Matrix<Double>& elemMat, EntityIterator& ent1, EntityIterator& ent2 );

  protected:
    Matrix<Double>& stabilParams_;
};

class FluidMechPlaneStiffNewtonInt_PV : public FluidMechInt
{
  public:
    FluidMechPlaneStiffNewtonInt_PV(Double density, \
        Double kinematicViscosity,Matrix<Double> & stabilParams, \
        bool movingMesh, std::string stabilType);
    virtual ~FluidMechPlaneStiffNewtonInt_PV();
    void CalcElementMatrix( Matrix<Double>& elemMat, EntityIterator& ent1, EntityIterator& ent2 );

  protected:
    Matrix<Double>& stabilParams_;
};

class FluidMechPlaneStiffNewtonInt_UQ : public FluidMechInt
{
  public:
    FluidMechPlaneStiffNewtonInt_UQ(Double density, \
        Double kinematicViscosity, Matrix<Double> & stabilParams, \
        bool movingMesh, std::string stabilType);
    virtual ~FluidMechPlaneStiffNewtonInt_UQ();
    void CalcElementMatrix( Matrix<Double>& elemMat, EntityIterator& ent1, EntityIterator& ent2 );

  protected:    
    Matrix<Double> & stabilParams_;
};

class FluidMechPlaneStiffNewtonInt_PQ : public FluidMechInt
{
  public:
    FluidMechPlaneStiffNewtonInt_PQ(Double density, \
        Double kinematicViscosity, Matrix<Double> & stabilParams, \
        bool movingMesh, std::string stabilType);
    virtual ~FluidMechPlaneStiffNewtonInt_PQ();
    void CalcElementMatrix( Matrix<Double>& elemMat, EntityIterator& ent1, EntityIterator& ent2 );

  protected:    
    Matrix<Double> & stabilParams_;
};

class FluidMechPlaneIntNewton_RhsUV : public FluidMechInt
{
  public:
    FluidMechPlaneIntNewton_RhsUV(Double density, \
        Double kinematicViscosity,Matrix<Double> & stabilParams, \
        bool movingMesh, std::string stabilType);
    virtual ~FluidMechPlaneIntNewton_RhsUV();
    void CalcElementMatrix( Matrix<Double>& elemMat, EntityIterator& ent1, EntityIterator& ent2 );

  protected:
    Matrix<Double>& stabilParams_;
};

class FluidMechPlaneIntNewton_RhsUQ : public FluidMechInt
{
  public:
    FluidMechPlaneIntNewton_RhsUQ(Double density, \
        Double kinematicViscosity,Matrix<Double> & stabilParams, \
        bool movingMesh, std::string stabilType);
    virtual ~FluidMechPlaneIntNewton_RhsUQ();
    void CalcElementMatrix( Matrix<Double>& elemMat, EntityIterator& ent1, EntityIterator& ent2 );

  protected:
    Matrix<Double>& stabilParams_;
};

class FluidMechAbsorbingFlow : public LinearSurfForm 
{
public:
  FluidMechAbsorbingFlow(Double kinematicViscosity ); 

  virtual ~FluidMechAbsorbingFlow();

    /// Calculation of RHS vector for double entries, i.e. transient and static 
    virtual void CalcElemVector( Vector<Double> & result,
                                 EntityIterator& ent );

protected:
  Double kinematicViscosity_;
  Matrix<Double> elemResult_;
};



} // namespace
#endif //FILE_FLUIDMECHSTIFFINT

