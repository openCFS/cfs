#ifndef FILE_CFS_ANSATZFCT_HH
#define FILE_CFS_ANSATZFCT_HH


namespace CoupledField {


  //! Base class representing an ansatz function space
  class AnsatzFct {

  public:
    typedef enum {CONST, LAGRANGE, LEGENDRE} AnsatzFctType;

    //! Constuctor
    AnsatzFct();

    //! Return type of ansatz functions
    AnsatzFctType GetType() { return type_; }

  protected:
    
    //! Type of ansatz functions used
    AnsatzFctType type_;

  };

  //! Comparison operator for AnsatzFct
  bool operator==( const AnsatzFct& a, const AnsatzFct& b );

  //! Negative comparison operator for AnsatzFct
  bool operator!=( const AnsatzFct& a, const AnsatzFct& b );

  //! Class for constant ansatz functions
  class ConstFct : public AnsatzFct {

  public:
    
    //! Constructor
    ConstFct();

  };

  //! Class for standard 1st and 2nd orer Lagrange functions
  class LagrangeFct : public AnsatzFct {

  public:

    //! Constructor
    LagrangeFct();
    
  };

  //! Class for hierarchic Legendre ansatz functions
  class LegendreFct { 

  public:

    //! Constructor
    LegendreFct();
    
  };

}



#endif
