#ifndef FILE_FLUIDMECHSHEARSTRESS
#define FILE_FLUIDMECHSHEARSTRESS

#include "Forms/linElastInt.hh"
#include "General/environment.hh"

#include "General/defs.hh"
#include "MatVec/matrix.hh"

namespace CoupledField {
class BaseMaterial;
class EntityIterator;
template <class TYPE> class Vector;
}  // namespace CoupledField

namespace CoupledField
{
  
  //! class for calculation of shear stresses
  template <class TYPE>
  class FluidMechShearStress : public linElastInt
  {
  public:

    /// Constructor
    FluidMechShearStress(BaseMaterial* matData, SubTensorType type);
  
    /// Destructor
    ~FluidMechShearStress();  
  
    /// calculates shear-stresses (vector notation)
    void CalcShearStressVec( Matrix<TYPE>& stressVec, \
                             EntityIterator& volEnt, \
                             const EntityIterator& surfEnt, \
                             const Vector<Double> outNormal);
    void CalcPressureForce( Matrix<TYPE>& pressureForce, \
                            EntityIterator& volEnt, \
                            const EntityIterator& surfEnt,
                            const Vector<Double> outNormal);

    /// in stress calculations, the actual velocities of the element is needed
    /*!
      \param velo (input) Matrix with velocities v of all nodes of actual element
      \f[ \left( \begin{array}{ccc} 
      v_{x1} &  v_{x2} &  v_{x3} \\
      v_{y1} &  v_{y2} &  v_{y3} \\
      \end{array}\right) \f]         
    */
    void SetActElemSol(Matrix<TYPE>& velo) {
      elemVelo_ = velo;};
    void SetActElemPresSol(Vector<TYPE>& presSol) {
      elemPres_ = presSol;};
    
  protected:  
  
    /// returns B 
    void CalcBMat(Matrix<Double> & bMat, UInt ip, const Matrix<Double> & ptCoord);

    /// velocities of all nodes of actual element
    Matrix<TYPE> elemVelo_;
    Vector<TYPE> elemPres_;
    Double density_, dynamicViscosity_, kinematicViscosity_;

  };
  

} //end namespace

#endif // FILE_FLUIDMECHSHEARSTRESS


