// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_LINSTRAINOP_2005
#define FILE_LINSTRAINOP_2005

#include "Forms/baseoperator.hh"
#include "Utils/tools.hh"
#include "Utils/StdVector.hh"

namespace CoupledField
{

  class Grid;
  struct Elem;
  class StdPDE;
  template<class TYPE> class NodeStoreSol;
  template<class TYPE> class Vector;
  template<class TYPE> class Matrix;
  
  //! Operator for calculating linear Strain fields
  
  //! This operator class calculates linear Strains, i.e. the 
  // differential operator B applied to the mechanical displacement u

  template<class TYPE>
  class LinStrainOp : public BaseOperator
  {

  
  public:

    //! Constructor
    /*!
      \param ptGrid (input) Pointer to grid
      \param ptPDE (input) Pointer to PDE
      \param ptEQN (input) Pointer to EQN
      \param potential (input) NodeStoreSol containing nodal potential
      \param solType (input) SolutionType of the potentialField
      \param isaxi (input) Flag for axi-symmetric geomtetry
    */
    LinStrainOp(Grid * ptGrid,
                StdPDE * ptPDE,
                shared_ptr<EqnMap> eqnMap,
                NodeStoreSol<TYPE> & displacement_,
                const SolutionType solType,
                bool isaxi=false);
 
    //! Destructor
    virtual ~LinStrainOp();
  
    //! Calculate element linear Strain
    /*!
      \param linStrain (output) Element vector of linear Strain field
      \param ptElem (input) Pointer to element
      \param LCoord (input) Local coordinates of evaluation point
    */
    virtual void CalcElemLinearStrain(SingleVector & linStrain,
                                   const Elem * ptElement,
                                   Matrix<Double> & LCoord);
  

    //! Calculate linear Strain for list of subdomains
    /*!
      \param linStrain (output) Vector 
      \param SD (input) Name of the subdomain
      \param LCoord (input) Local coordinates of evalutation point
    */
    virtual void CalcSDLinearStrain(SingleVector & linStrain,
                                 const StdVector<RegionIdType> & SD,
                                 Matrix<Double> & lCoord);

    void calcBMat(Matrix<Double> & bMat,BaseFE * ptelem, UInt ip, 
                  Matrix<Double> & ptCoord);
                                                       

  protected:
  
    //! StoreSolution containing potentialfield
    NodeStoreSol<TYPE> * displacement_;

    //! Soltution type of the potential field
    SolutionType solType_;

    StdPDE * ptPDE_;
  };


} // end of namespace

#endif
