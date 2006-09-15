#ifndef FILE_BUBBLESTIFF_INT
#define FILE_BUBBLESTIFF_INT

#include "baseForm.hh"

namespace CoupledField
{

  //! Class for calculation  element mass matrix
  class BubbleStiffInt : public BaseForm
  {
  public:

    // Constructor
    BubbleStiffInt(const Double aDensity,
		   Double densityforbubble,
		   Double sonicVel, 
		   Double viscosity,
		   Double bubbleDensity,
		   Double frequency,
		   bool axi=false );

    // Destructor
    virtual ~BubbleStiffInt();

    // Calculation of element matrix
    void CalcElementMatrix( Matrix<Double>& elemMat,
                            EntityIterator& ent1, 
                            EntityIterator& ent2 );
      
    //! set additional multiplicative factor of mass matrix
    void SetFactor(Double aFactor) 
    {factor_ = aFactor;};
   
    //! Set element Numbers and radius values
    void SetValues( StdVector<UInt> elemNumbers, Vector<Double> * radius,
		    Vector<Double> * radiusDeriv );
 

  protected:
    

  
  private:

    Double density_;          //!< multiplicative value for mass integrator
    Double factor_;           //!< yet another multiplicative value for mass integrator
 
    //! vector for element number
    StdVector<UInt> elemNumbers_;

    //! pointer to vector with element radius values
    Vector<Double> * radius_;

    //! pointer to vector with element radius values derivatives
    Vector<Double> * radiusDeriv_;
    
    //! map from global element number to index position in elemNumbers/radius
    std::map<UInt, UInt> indexMap_;

    //! parameters to compute bubble-dependend coefficients
    Double sonicVel_, densityforbubble_, viscosity_, bubbleDensity_;

    //! parameters to determine when bubbleterms should be added
    Double frequency_;
  
  };





}

#endif // FILE_MASSINT
