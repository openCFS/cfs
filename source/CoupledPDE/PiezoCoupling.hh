// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_PIEZOCOUPLING_HH
#define FILE_PIEZOCOUPLING_HH

#include "BasePairCoupling.hh"
#include "Utils/elemstoresol.hh"
#include "Utils/nodestoresol.hh"
#include "PDE/SinglePDE.hh"



namespace CoupledField
{

  // Forward declarations
  class BaseForm;
  class BaseMaterial;
  class SinglePDE;

  //! Implements the definition of pairwise piezo-coupling
  class PiezoCoupling : public BasePairCoupling
  {
  public:
    //! Constructor
    //! \param pde1 pointer to first coupling PDE
    //! \param pde2 pointer to second coupling PDE
    //! \param paramNode pointer to "couplinglist/direct/piezoDirect" element
    PiezoCoupling( SinglePDE *pde1, SinglePDE *pde2, ParamNode * paramNode );

    //! Destructor
    virtual ~PiezoCoupling();

    //! Trigger calculation of postprocessing results
    void CalcResults( shared_ptr<BaseResult> result );

//     template <class TYPE>
//     void calcMaterialMatrices(Matrix<TYPE> &sMat, 
//                               Matrix<TYPE>&cMat,
//                               Matrix<TYPE> &pMat,
//                               Matrix<Double> *matDat,
//                               Matrix<Complex> *complexMatDat);

 

    //! Gathers all information concerning nonlinear computations
    void ReadPiezoNonLin();

    //!
    void GetNonlinMaterialTensor( Matrix<Double>& matTensor, 
                                  Vector<Double>& elecD,
                                  std::string matTensorType,
                                  BaseMaterial* matMech,
                                  BaseMaterial* matElec,
                                  BaseMaterial* matCouple,
                                  SubTensorType subTensorType,
                                  Vector<Double>& elecField,
                                  Vector<Double>& mechStrain,
                                  Vector<Double>& elecFieldPrev,
                                  Vector<Double>& mechStrainPrev,
                                  EntityIterator& ent );
    
  protected:

    //! Definition of the (bi)linear forms
    void DefineIntegrators();
   
    //! Define available results
    void DefineAvailResults();

 
    // Data section
    bool hasOutput_;

    //! flag, which is true if piezo coupling is nonlinear, i.e. 
    //! if we have a nonlinear piezoelectric coupling form
    //! it is not the same as nonLin_ which is already set true if one of the
    //! singlePDEs is nonlinear
    bool nonLinPiezoCoupling_;

  private:

    //! compute normalized irreversible strain
    void ComputeSirr( Vector<Double>& VecSirr, SubTensorType type,
                      UInt dirP, Double ctP, 
                      BaseMaterial* matMech );

    //!
    void ComputeDiffCouplingTensor( Matrix<Double>& dMat, 
                                    Vector<Double>& actE,
                                    Vector<Double>& prevE,
                                    Vector<Double>& actSirr,
                                    Vector<Double>& prevSirr,
                                    Directions dirP,
                                    SubTensorType subTensorType );
    // Postprocession section

    //! computes stresses, strain, i.e. \sigma = cBu + e \grad \phi
    template <class TYPE>
    void CalcStressStrain( shared_ptr<BaseResult> result );

    //! computes irreversibel strain
    void CalcStrainIrr( shared_ptr<BaseResult> result );

    //! calculate Charges
    template <class TYPE>
    void CalcCharges( shared_ptr<BaseResult> result );

    //! calculate DField
    template <class TYPE>
    void CalcDField( shared_ptr<BaseResult> result );

    // Data section

    //! flag indicating use of complex material parameters
    bool hasComplexMatParams_;

  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class PiezoCoupling
  //! 
  //! \purpose 
  //! 
  //! \collab 
  //! 
  //! \implement 
  //! 
  //! \status In use
  //! 
  //! \unused 
  //! 
  //! \improve
  //! 

#endif

} // end of namespace

#endif
