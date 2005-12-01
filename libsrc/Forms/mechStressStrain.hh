#ifndef FILE_MECHSTRESSSTRAIN_04
#define FILE_MECHSTRESSSTRAIN_04

#include <Elements/basefe.hh>
#include <DataInOut/MaterialData.hh>

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
    MechStressStrain(BaseFE * aptelem, MaterialData & matData);

    /// Constructor
    MechStressStrain(MaterialData & matData);
  
    /// Destructor
    virtual ~MechStressStrain();  
  
    /// calculates Piola-Kirchoff-stresses (vector notation)
    void CalcStressVec(Vector<TYPE>& stressVec, UInt ip, 
                       Matrix<Double> & ptCoord);  

    // calculate linear part of Green Lagrangian Strain tensor
    void CalcStrainVec(Vector<TYPE>& strainVec, UInt ip, 
                       Matrix<Double> & ptCoord);  
    
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

    /// returns nr. of degrees of freedom
    virtual UInt getNrDofs()=0;  

    /// displacement of all nodes of actual element
    Matrix<TYPE> elemDisp_;

  };
  

  //! 3D Cauchy stresses, linear strains
  template <class TYPE>
  class MechStressStrain3D : public MechStressStrain<TYPE>
  {
  public:

    /// Constructor
    MechStressStrain3D(BaseFE * aptelem, MaterialData & matData);

    /// Constructor
    MechStressStrain3D(MaterialData & matData);
  
    /// Destructor
    virtual ~MechStressStrain3D();  
  
  protected:  
    /// returns D - matrix (material matrix)
    virtual void calcDMat(Matrix<Double> & dMat);

    /// returns nr. of degrees of freedom
    virtual UInt getNrDofs(){return 3;};

    /// returns the size of the material d-matrix
    virtual UInt getDimD(){return 6;};

  };
  



  // second part: regarding internal stresses
  template <class TYPE>
  class MechStressStrainPlaneStrain : public MechStressStrain<TYPE>
  {
  public:
    /// Constructor
    MechStressStrainPlaneStrain(BaseFE * aptelem, MaterialData & matData);

    /// Constructor
    MechStressStrainPlaneStrain(MaterialData & matData);
  
    /// Destructor
    virtual ~MechStressStrainPlaneStrain();  
  
  protected:  
    /// returns D - matrix (material matrix)
    virtual void calcDMat(Matrix<Double> & dMat);

    /// returns nr. of degrees of freedom
    virtual UInt getNrDofs(){return 2;};  

    /// returns the size of the material d-matrix
    virtual UInt getDimD(){return 3;};

  };
  

  // second part: regarding internal stresses
  template <class TYPE>
  class MechStressStrainAxi : public MechStressStrain<TYPE>
  {
  public:
    /// Constructor
    MechStressStrainAxi(BaseFE * aptelem, MaterialData & matData);

    /// Constructor
    MechStressStrainAxi(MaterialData & matData);
  
    /// Destructor
    virtual ~MechStressStrainAxi();  
  
  protected:  
    /// returns D - matrix (material matrix)
    virtual void calcDMat(Matrix<Double> & dMat);

    /// returns nr. of degrees of freedom
    virtual UInt getNrDofs(){return 2;};  

    /// returns the size of the material d-matrix
    virtual UInt getDimD(){return 4;};
  };
  

  // Explicite template instantiation
  template class MechStressStrain<Double>;
  template class MechStressStrain<Complex>;

  template class MechStressStrain3D<Double>;
  template class MechStressStrain3D<Complex>;

  template class MechStressStrainPlaneStrain<Double>;
  template class MechStressStrainPlaneStrain<Complex>;

  template class MechStressStrainAxi<Double>;
  template class MechStressStrainAxi<Complex>;
} //end namespace

#endif // FILE_XXX


