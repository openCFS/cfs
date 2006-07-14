#ifndef FILE_FLATSHELLSTIFFINT
#define FILE_FLATSHELLSTIFFINT

#include "FlatShellInt.hh"

namespace CoupledField
{
     //! Class for flat shell Integrator
    class FlatShellStiffInt : public FlatShellInt
 {
    public:
    
    //! Constructor with pointer to normal material
    FlatShellStiffInt(BaseMaterial *matData);
    
    //! Constructor with pointer to composite material
    FlatShellStiffInt( Composite *composite);
    
    //!  Destructor
    virtual ~FlatShellStiffInt();

    //!returns dimension of D matrix of bending stiffness
    virtual Integer getDimD(){return 8;};
  
    //!returns nr. of degrees of freedom
    virtual Integer getNrDofs(){return 6;};
   				    
   //! Enumeration for differentiation of the three mechanical parts Membrane, Bending, Shear
   typedef enum {MEMBRANE,BENDING,SHEAR,COUPLED1,COUPLED2,PIEZO} SubPartType;
   
   //! Function for calculation bdb matrix 
   void CalcElementMatrix( Matrix<Double>& elemMat,
			   EntityIterator& ent1, 
			   EntityIterator& ent2 );
   
   //! returns B - matrix for BDB
   void calcBMat(Matrix<Double> & bMat, Integer ip, 
                 BaseFE* elem, Matrix<Double> & ptCoord, SubPartType parts );
    
    StdVector<SubPartType> parts_;
    StdVector<bool> reducedPart_;
    
  protected:
    
    //!calculates the material Matrix for Normal Material
    void calcDMat( Matrix<Double> &dMat);
    
    //!calculates the material Matrix for Composite Material
    void calcDMatComposite( Matrix<Double> &dMat );
                   
    //Parametres of the lamina structures
    Vector<Double> z_;
    Vector<Double> orAngle_;
       
 };

}//end of namespace

#endif 
