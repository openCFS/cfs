// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_GRADIENTFIELDOP_2004
#define FILE_GRADIENTFIELDOP_2004

#include "Forms/baseoperator.hh"
#include "Utils/tools.hh"
#include "Utils/StdVector.hh"

namespace CoupledField
{

  class Grid;
  struct Elem;
  template<class TYPE> class NodeStoreSol;
  template<class TYPE> class Vector;
  template<class TYPE> class Matrix;
  
  //! Operator for calculating gradient fields
  
  //! This operator class calculates elemt gradient fields like electric
  //! field or electric flux density at arbitrary points for a given
  //! potential field

  template<class TYPE>
  class GradientFieldOp : public BaseOperator
  {

  
  public:

    //! Constructor
    /*!
      \param ptGrid (input) Pointer to grid
      \param ptPDE (input) Pointer to PDE
      \param ptEQN (input) Pointer to EQN
      \param potential (input) NodeStoreSol containing nodal potential.
      \param solType (input) SolutionType of the potentialField
      \param isaxi (input) Flag for axi-symmetric geomtetry
      \param coordUpdate (input) flag indicating updated Lagrangian 
                                 formulation
    */
    GradientFieldOp(Grid * ptGrid,
                    StdPDE * ptPDE,
                    shared_ptr<EqnMap> eqnMap,
                    NodeStoreSol<TYPE> & potential,
                    shared_ptr<AnsatzFct> fctType,
                    bool isaxi=false,
                    bool coordUpdate = false );

    /** Alternative constructor where the element solution is given
     * in CalcElemGradField() as parameter */
    GradientFieldOp(Grid * ptGrid,
                    StdPDE * ptPDE,
                    shared_ptr<AnsatzFct> fctType,
                    bool isaxi=false,
                    bool coordUpdate = false );


    //! Destructor
    virtual ~GradientFieldOp();
  
    //! Calculate element gradient field
    /*!
      \param elemField (output) Element vector of gradient field
      \param ent (input) EntityIterator pointing to current element
      \param lCoord (input) Local coordinates of evaluation point
      \param factor (input) Scaling factor (e.g. permittivity for E-Field)
      \param dof if the element describes a vector field, dof (1-based) selects the scalar field
      \param elem_data optionally we give the read data and ignore potential
    */
    void CalcElemGradField(Vector<TYPE> & elemField,
                                   const EntityIterator& ent,
                                   const Vector<Double> & lCoord,
                                   const Double factor,
                                   const UInt dof = 1,
                                   const SingleVector* elem_data = NULL);
  

    //! Calculate electric field for list of subdomains
    /*!
      \param elemField (input) Vector containing
      \f[ elemFIeld[0] = \left( elemField_{1,x}, elemField_{2,x}, \cdots \right),
      elemField[1] = \left( elemField_{1,y}, elemField_{2,y}, \cdots \right), \cdots 
      \f]
      \param SD (input) Name of the subdomain
      \param LCoord (input) Local coordinates of evalutation point
    */
    void CalcSDGradField(Vector<TYPE> & elemField,
                                 const StdVector<RegionIdType> & SD,
                                 const Vector<Double> & lCoord,
                                 const Vector<Double> & factors);
                                                       

  protected:
  
    //! StoreSolution containing potentialfield
    NodeStoreSol<TYPE> * potential_;

    /** the ansatz functions from the ResultInfo */
    shared_ptr<AnsatzFct> fctType_;
  };

} // end of namespace

#endif
