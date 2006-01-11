#ifndef FILE_ACOUENERGYOP
#define FILE_ACOUENERGYOP

#include "Forms/baseoperator.hh"



namespace CoupledField {

  // Forward declaration of classes
  class Grid;
  class Elem;
  template<class TYPE> class Vector;
  template<class TYPE> class Matrix;
  
  
  //! This operator class calculates the acoustic power density per
  //! element for a acoustic PDE formulation in velocity potential

  class AcouPowerDensityOp : public BaseOperator
  {

  
  public:

    //! Constructor
    /*!
      \param ptGrid (input) Pointer to grid
      \param ptPDE  (input) Pointer to PDE
      \param ptEQN  (input) Pointer to EQN
      \param isaxi  (input) Flag for axi-symmetric geomtetry
    */
    AcouPowerDensityOp(Grid * ptGrid,
                       StdPDE * ptPDE,
                       NodeEQN * ptEQN,
                       Boolean isaxi=FALSE);
    
    //! Destructor
    virtual ~AcouPowerDensityOp();
  
    //! Calculate acoustic power density for element
    /*!
      \param elemPD (output) Element vector of acoustic power density
      \param ptElem (input) Pointer to element
      \param density (input) multiplicative factor
    */
    virtual void CalcElemPD(Vector<Double> & elemPD,
                                const Elem * ptElement,
                                const Double density);
    

    //! Calculate acoustic power density for elements in subdomains
    /*!
      \param elemPD (output) Vector containing
      \param SD (input) Identifier of the subdomain
      \param density (input) multiplicative factor
    */
    virtual void CalcSubdomPD(Vector<Double> & elemPD,
                              const StdVector<RegionIdType> & SD,
                              const Double density);
                                                       

  protected:
  
    //NodeStoreSol<Double> * sol_;
    //NodeStoreSol<Double> * solDeriv1_;
    //SolutionType solType_; //!< Soltution type

    Boolean isaxi_;
  };


} // end of namespace

#endif
