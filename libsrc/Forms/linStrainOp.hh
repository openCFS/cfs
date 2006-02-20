#ifndef FILE_LINSTRAINOP_2005
#define FILE_LINSTRAINOP_2005

#include "Forms/baseoperator.hh"
#include "Utils/tools.hh"
#include "Utils/StdVector.hh"

#include "olas.hh"

namespace CoupledField
{

  class Grid;
  class Elem;
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
                NodeEQN * ptEQN,
                NodeStoreSol<TYPE> & displacement_,
                const SolutionType solType,
                Boolean isaxi=FALSE);
 
    //! Destructor
    virtual ~LinStrainOp();
  
    //! Calculate element linear Strain
    /*!
      \param linStrain (output) Element vector of linear Strain field
      \param ptElem (input) Pointer to element
      \param LCoord (input) Local coordinates of evaluation point
    */
    virtual void CalcElemLinearStrain(CFSVector & linStrain,
                                   const Elem * ptElement,
                                   Matrix<Double> & LCoord);
  

    //! Calculate linear Strain for list of subdomains
    /*!
      \param linStrain (output) Vector 
      \param SD (input) Name of the subdomain
      \param LCoord (input) Local coordinates of evalutation point
    */
    virtual void CalcSDLinearStrain(CFSVector & linStrain,
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
