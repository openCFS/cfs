#ifndef FILE_FLUIDMECHSTIFFINT
#define FILE_FLUIDMECHSTIFFINT

#include "fluidMechInt.hh"

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

} // namespace
#endif //FILE_FLUIDMECHSTIFFINT

