#ifndef FILE_MECHSTRESSSTRAIN_04
#define FILE_MECHSTRESSSTRAIN_04

#include <Elements/basefe.hh>
#include <DataInOut/MaterialData.hh>

#include <Forms/linElastInt.hh>
#include <General/environment.hh>

namespace CoupledField
{
  
  //! class for calculation of mechanical stresses
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
    void CalcStressVec(Vector<Double>& StressVec, Integer ip, Matrix<Double> & ptCoord);  

    /// in stress calculations, the actual displacement of the element is needed
    /*!
      \param disp (input) Matrix with displacement d of all nodes of actual element
      \f[ \left( \begin{array}{ccc} 
      d_{x1} &  d_{x2} &  d_{x3} \\
      d_{y1} &  d_{y2} &  d_{y3} \\
      \end{array}\right) \f]         
    */
    virtual void SetActElemSol(Matrix<Double>& disp) {
      ENTER_FCN( "MechStressStrain::SetActElemSol" );
      elemDisp_ = disp;};

  protected:  
    /// returns D - matrix (material matrix)
    //  virtual void calcDMat(Matrix<Double> & dMat);
  
    /// returns B 
    virtual void calcBMat(Matrix<Double> & bMat, Integer ip, Matrix<Double> & ptCoord);

    /// returns nr. of degrees of freedom
    virtual Integer getNrDofs()=0;  

    /// displacement of all nodes of actual element
    Matrix<Double> elemDisp_;

  };
  

  //! 3D Cauchy stresses, linear strains
  class MechStressStrain3D : public MechStressStrain
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
    virtual Integer getNrDofs(){return 3;};

    /// returns the size of the material d-matrix
    virtual Integer getDimD(){return 6;};

  };
  



  // second part: regarding internal stresses
  class MechStressStrainPlaneStrain : public MechStressStrain
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
    virtual Integer getNrDofs(){return 2;};  

    /// returns the size of the material d-matrix
    virtual Integer getDimD(){return 3;};

  };
  

  // second part: regarding internal stresses
  class MechStressStrainAxi : public MechStressStrain
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
    virtual Integer getNrDofs(){return 2;};  

    /// returns the size of the material d-matrix
    virtual Integer getDimD(){return 4;};
  };
  

} //end namespace

#endif // FILE_XXX


