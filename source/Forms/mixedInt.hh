#ifndef FILE_MIXEDINT_1
#define FILE_MIXEDINT_1

#include "baseForm.hh"
#include "massInt.hh"

namespace CoupledField
{

  //! Class for calculation  element mass matrix
  class MassMixedInt_PP : public BaseForm
  {
  public:

    // Constructor
    MassMixedInt_PP(const Double aDactor, bool axi=false, bool coordUpdate = false );

    // Destructor
    virtual ~MassMixedInt_PP();

    // Calculation of element matrix
    void CalcElementMatrix( Matrix<Double>& elemMat,
                            EntityIterator& ent1, 
                            EntityIterator& ent2 );
      

  protected:
  
  private:

    Double factor_;          //!< multiplicative value for mass integrator
    UInt nrDofsPerNode_;   //!< degrees of freedom per node
  };

  //! Class for calculation  element mass matrix
  class MassMixedInt_VV : public BaseForm
  {
  public:

    // Constructor
    MassMixedInt_VV(const Double aDactor, bool axi=false, bool coordUpdate = false );

    // Destructor
    virtual ~MassMixedInt_VV();

    // Calculation of element matrix
    void CalcElementMatrix( Matrix<Double>& elemMat,
                            EntityIterator& ent1, 
                            EntityIterator& ent2 );
      

  protected:
  
  private:

    Double factor_;          //!< multiplicative value for mass integrator
    UInt nrDofsPerNode_;   //!< degrees of freedom per node
  };



  //! Class for calculation  element stiffness matrix 
  class StiffMixedInt_KPV : public BaseForm
  {
  public:
    
    /// Constructor
    StiffMixedInt_KPV(Double laplVal, bool axi=false, 
		    bool coordUpdate = false );
    
    /// 
    virtual ~StiffMixedInt_KPV();
    
    /// Calculation of stiffmess matrix
    void CalcElementMatrix( Matrix<Double>& elemMat,
			    EntityIterator& ent1, 
			    EntityIterator& ent2 );
    
    
    
  private:
    /// multiplicative value for laplace integration 
    Double factor_;
  };




  //! Class for calculation  element stiffness matrix 
  class StiffMixedInt_KVP : public BaseForm
  {
  public:
    
    /// Constructor
    StiffMixedInt_KVP(Double laplVal, bool axi=false, 
		    bool coordUpdate = false );
    
    /// 
    virtual ~StiffMixedInt_KVP();
    
    /// Calculation of stiffmess matrix
    void CalcElementMatrix( Matrix<Double>& elemMat,
			    EntityIterator& ent1, 
			    EntityIterator& ent2 );
    
    
    
  private:
    /// multiplicative value for laplace integration 
    Double factor_;
  };




}

#endif // FILE_MIXEDINT
