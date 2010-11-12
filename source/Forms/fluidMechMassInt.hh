#ifndef FILE_FLUIDMECHMASSFULLNEWTONINT
#define FILE_FLUIDMECHMASSFULLNEWTONINT

#include "fluidMechInt.hh"

namespace CoupledField
{

//**********************************************************
//************** Newton Stabilized FEM ****************
//**********************************************************

class FluidMechPlaneMassNewtonInt_UV : public FluidMechInt
{
  public:
    FluidMechPlaneMassNewtonInt_UV(Double density,
        Double kinematicViscosity,Matrix<Double> & stabilParams,
        bool movingMesh, std::string stabilType);
    virtual ~FluidMechPlaneMassNewtonInt_UV();
    void CalcElementMatrix( Matrix<Double>& elemMat, EntityIterator& ent1, EntityIterator& ent2 );

  protected:
    Matrix<Double> & stabilParams_;
};

class FluidMechPlaneMassNewtonInt_UQ : public FluidMechInt
{
  public:
    FluidMechPlaneMassNewtonInt_UQ(Double density,
        Double kinematicViscosity,Matrix<Double> & stabilParams,
        bool movingMesh, std::string stabilType);
    virtual ~FluidMechPlaneMassNewtonInt_UQ();
    void CalcElementMatrix( Matrix<Double>& elemMat, EntityIterator& ent1, EntityIterator& ent2 );

  protected:
    Matrix<Double> & stabilParams_;
};

} // namespace
#endif // FILE_FLUIDMECHMASSFULLNEWTONINT
