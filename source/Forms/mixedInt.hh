#ifndef FILE_MIXEDINT_1
#define FILE_MIXEDINT_1

#include "baseForm.hh"
#include "massInt.hh"
#include "Forms/linSurfForm.hh"

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


  //! Class for calculation  of ABC in mixed formulation
  class ABC_MixedInt : public BaseForm
  {
  public:

    // Constructor
    ABC_MixedInt(const Double aDactor, bool axi=false, bool coordUpdate = false );

    // Destructor
    virtual ~ABC_MixedInt();

    // Calculation of element matrix
    void CalcElementMatrix( Matrix<Double>& elemMat,
                            EntityIterator& ent1, 
                            EntityIterator& ent2 );
      

  protected:
  
  private:

    Double factor_;          //!< multiplicative value for mass integrator
    UInt nrDofsPerNode_;   //!< degrees of freedom per node
  };


  //! Class implementing an inhomogeneous Neumann integrator for acoustic field
  class LinSurfVelocity : public LinearSurfForm {

  public:
    
    //! Standard constructor
    LinSurfVelocity( std::string amplitudeStr, std::string phaseStr,
		     Double val, bool isaxi );
    
    //! Destructor
    ~LinSurfVelocity();

    /// Calculation of RHS vector for double entries, i.e. transient and static 
    void CalcElemVector( Vector<Double> & elemVec,
                         EntityIterator& ent );

    /// Calculation of RHS vector for complex entries, i.e. harmonic
    void CalcElemVector( Vector<Complex> & elemVec,
                         EntityIterator& ent );
    
    void SetAmplitude(const std::string& value) { amplitude_ = value; }
    
  protected:

    /// comprises part of CalcElemVector, which is equal for double and complex
    void PrepareElemVector( Vector<Double> & elemVec,
                            EntityIterator& ent );
    
    //! amplitude string
    std::string amplitude_;
    
    //! phase string
    std::string phase_;

    //! material factor
    Double matFactor_;
  };


}

#endif // FILE_MIXEDINT
