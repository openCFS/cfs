#ifndef FILE_NLINMASSINT_1
#define FILE_NLINMASSINT_1

#include "baseForm.hh"

namespace CoupledField
{

  //! Class for calculation  element mass matrix
  class NlinMassInt : public BaseForm
  {
  public:

    // Constructor
    NlinMassInt(const Double aDensity, const UInt nrDofsPerNode=1, 
            bool axi=false, bool coordUpdate = false );

    // Destructor
    virtual ~NlinMassInt();

    // Calculation of element matrix
    void CalcElementMatrix( Matrix<Double>& elemMat,
                            EntityIterator& ent1, 
                            EntityIterator& ent2 );
      

    //! for fractional damping model, is called in AcousticPDE::DefineIntegrators
    void SetFracDamping() 
    {isFracDamping_ = true;};

    //! for fractional damping model, is called in AcousticPDE::DefineIntegrators
    void UnsetFracDamping() 
    {isFracDamping_ = false;};
  
    //! set additional multiplicative factor of mass matrix
    void SetFactor(Double aFactor) 
    {factor_ = aFactor;};
   
  protected:
    
  
  private:

    Double density_;          //!< multiplicative value for mass integrator
    Double factor_;           //!< yet another multiplicative value for mass integrator
    UInt nrDofsPerNode_;   //!< degrees of freedom per node
  
  };





}

#endif // FILE_MASSINT
