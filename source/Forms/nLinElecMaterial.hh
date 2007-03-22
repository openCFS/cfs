// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_NLINELECMATERIAL_06
#define FILE_NLINELECMATERIAL_06

#include <Elements/basefe.hh>
#include <Materials/baseMaterial.hh>
#include <Forms/gradfieldop.hh>
#include <General/environment.hh>
#include <Utils/ApproxData.hh>
#include <Forms/bdbInt.hh>


namespace CoupledField
{

  // =============================================================================
  // 3D nonlinear material mechanics
  // =============================================================================
  
  
  // class for calculation of 3d nonlinear elasticity
  // first part: nonlinear B-matrix

  class nLinElec3dInt_Material : public BDBInt {
    
  public:
      
    /// Constructor
    nLinElec3dInt_Material(ApproxData *nlinFnc,  
                           BaseMaterial* matData, 
                           SubTensorType type = FULL);
    
    /// Destructor
    ~nLinElec3dInt_Material();  

    //! set objects for computation of E-field
    void Set4NonLinMaterial(Grid* ptGrid, 
                            StdPDE* ptPDE,
                            shared_ptr<EqnMap> eqnMap,
                            shared_ptr<ResultInfo> result);
  
  protected:  
    /// returns D - matrix for BDB
    virtual void calcDMat(Matrix<Double> & dMat, UInt ip, 
                          Matrix<Double> & ptCoord);
    
    /// returns B - matrix for BDB
    virtual void calcBMat(Matrix<Double> & bMat, UInt ip, 
                          Matrix<Double> & ptCoord );
    
    /// returns dimension of D matrix
    virtual UInt getDimD(){return dim_;};
    
    /// returns nr. of degrees of freedom
    virtual UInt getNrDofs(){return 1;};
    
    void CalcElementMatrix( Matrix<Double>& elemMat,
                            EntityIterator& ent1, 
                            EntityIterator& ent2 );

    //! Query material type for \f$D\f$ tensor
    virtual MaterialType getDMaterialType() { return ELEC_PERMITTIVITY;};

   
  private:
    ApproxData *nLinFnc_;

    /// scalar electric potential of all nodes of actual element
    Vector<Double> elemPot_;

    GradientFieldOp<Double> * EfieldOp_;

    // space dimension
    UInt dim_;

    // local copy of entity iterator
    EntityIterator ent1_;
      
    
  };
}


#endif // FILE_XXX
