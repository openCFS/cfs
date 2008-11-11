#ifndef FILE_FLUIDMECHSHEARSTRESS
#define FILE_FLUIDMECHSHEARSTRESS

#include <Elements/basefe.hh>
#include <Materials/baseMaterial.hh>

#include <Forms/linElastInt.hh>
#include <General/environment.hh>

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
    void CalcShearStressVec(Vector<TYPE>& stressVec, UInt ip, 
                            EntityIterator& ent);  

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
    
  protected:  
  
    /// returns B 
    void calcBMat(Matrix<Double> & bMat, UInt ip, 
                          Matrix<Double> & ptCoord);

    /// velocities of all nodes of actual element
    Matrix<TYPE> elemVelo_;
    Double density_, dynamicViscosity_, kinematicViscosity_;

  };
  

} //end namespace

#endif // FILE_FLUIDMECHSHEARSTRESS


