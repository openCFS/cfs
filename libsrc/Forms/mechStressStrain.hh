#ifndef FILE_MECHSTRESSSTRAIN_04
#define FILE_MECHSTRESSSTRAIN_04

#include <Elements/basefe.hh>
#include <Materials/baseMaterial.hh>

#include <Forms/linElastInt.hh>
#include <General/environment.hh>

namespace CoupledField
{
  
  //! class for calculation of mechanical stresses
  template <class TYPE>
  class MechStressStrain : public linElastInt
  {
  public:

    /// Constructor
    MechStressStrain(BaseMaterial* matData, SubTensorType type);
  
    /// Destructor
    virtual ~MechStressStrain();  
  
    /// calculates Piola-Kirchoff-stresses (vector notation)
    void CalcStressVec(Vector<TYPE>& stressVec, UInt ip, 
                       EntityIterator& ent);  

    // calculate linear part of Green Lagrangian Strain tensor
    void CalcStrainVec(Vector<TYPE>& strainVec, UInt ip, 
                       EntityIterator& ent);  
    
    /// in stress calculations, the actual displacement of the element is needed
    /*!
      \param disp (input) Matrix with displacement d of all nodes of actual element
      \f[ \left( \begin{array}{ccc} 
      d_{x1} &  d_{x2} &  d_{x3} \\
      d_{y1} &  d_{y2} &  d_{y3} \\
      \end{array}\right) \f]         
    */
    virtual void SetActElemSol(Matrix<TYPE>& disp) {
      ENTER_FCN( "MechStressStrain::SetActElemSol" );
      elemDisp_ = disp;};

  protected:  
    /// returns D - matrix (material matrix)
    //  virtual void calcDMat(Matrix<Double> & dMat);
  
    /// returns B 
    virtual void calcBMat(Matrix<Double> & bMat, UInt ip, 
                          Matrix<Double> & ptCoord);

    /// displacement of all nodes of actual element
    Matrix<TYPE> elemDisp_;

  };
  

} //end namespace

#endif // FILE_XXX


