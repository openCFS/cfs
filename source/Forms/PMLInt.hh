// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_PMLINT
#define FILE_PMLINT

#include "baseForm.hh"
#include "pmlBasics.hh"
#include "linElastInt.hh"

namespace CoupledField
{


  /// Class for calculation  element stiffness/mass matrix for PML formulation
  
  class PMLInt : public BaseForm
  {
  public:
    
    //! Constructor
    PMLInt(std::string type, Double factor, std::string dampingTypePML, 
           Double damp, bool axi=false);
    
    /// 
    virtual ~PMLInt();
    
    //! Calculation of stiffmess matrix
    void CalcElementMatrix( Matrix<Complex>& elemMat,
                            EntityIterator& ent1, 
                            EntityIterator& ent2 );
    
    //! set min/max of x,y,z coordinates form where PML starts and ends
    void SetPosPML(Matrix<Double> & inner, Matrix<Double> & outer);
    
    //! set number of dofs per node
    void SetNrDofs( UInt numDofs ) {
      nrDofsPerNode_ = numDofs;
    };

    //! set actual solution
    void SetActElemSol(Matrix<Double>& disp){};
    
    
  private:
    
    //! Calculation of stiffmess matrix
    void CalcElementMatrixStiff(Matrix<Double> & ptCoord, Matrix<Complex> & elemMat);
    
    //! Calculation of mass matrix
    void CalcElementMatrixMass(Matrix<Double> & ptCoord, Matrix<Complex> & elemMat);
    
    //! object containing standard PML methods
    PMLBasics *pmlFnc_;
    
    //! multiplicative factor for forms
    Double formsFactor_;

    //! number of dofs per node
    UInt nrDofsPerNode_;
    
};


  /// Class for calculation mechanical element stiffness matrix for PML formulation
  
  class MechPMLInt : public linElastInt
  {
  public:
    
    //! Constructor
    MechPMLInt(std::string type,  BaseMaterial* matData,  
               std::string dampingTypePML, Double damp, 
               SubTensorType tensorType = FULL);
    
    /// 
    virtual ~MechPMLInt();
    
    //! Calculation of stiffmess matrix
    void CalcElementMatrix( Matrix<Complex>& elemMat,
                            EntityIterator& ent1, 
                            EntityIterator& ent2 );
    
    //! set min/max of x,y,z coordinates form where PML starts and ends
    void SetPosPML(Matrix<Double> & inner, Matrix<Double> & outer);
    
    void SetActElemSol(Matrix<Double>& disp){};
    
    
  private:
    
    //! Calculation of stiffmess matrix
    void calcBMatPML( Matrix<Double>& bMat, Vector<Double>& CoordAtIP,
                      Matrix<Complex> & bMatComplex, Complex jacDetC);
    
    //! object containing standard PML methods
    PMLBasics *pmlFnc_;
};

}

#endif // FILE_PMLINT
