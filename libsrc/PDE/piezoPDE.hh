#ifndef FILE_PIEZOPDE
#define FILE_PIEZOPDE


#include "basePDE.hh"

 
namespace CoupledField
{

  //! Class for piezoelectric simulations

  //! Detailed description to follow ...
  class PiezoPDE: public BasePDE
  {

  public:

    //! Constructor
    
    //! The constructor relies on the constructor of the BasePDE class and
    //! only preforms additional operations necessary for the special piezo
    //! case.
    //! The following default values are used: incStopCrit = 1e-2,
    //! residualStopCrit = 1e-3.
    //! \param aGrid pointer to grid object
    //! \param aBCs pointer to object describing boundary conditions
    //! \param aTimeFunc pointer to object of class TimeFunc
    //! \param aInFile pointer to object of class FileType; contains input
    //!        data.
    //! \param aOutFile pointer to object of class WriteResults; needed for
    //!        outputting data
    PiezoPDE( Grid *aGrid, BCs *aBCs, TimeFunc *aTimeFunc, FileType *aInFile,
	      WriteResults *aOutFile );

    //!  Destructor
    virtual ~PiezoPDE() {};


    //! define all (bilinearform) integrators needed for this pde
    virtual void DefineIntegrators(const Integer level);

    /// returns a stiffness integrator appropriate to the actual problem (e.g. 3D)
    BaseForm * GetStiffIntegrator(MaterialData& actSDMat, Boolean reducedInt=FALSE);

    //! prepare for correct time stepping
    /*! \param dt time step  */
    virtual void InitTimeStepping(const Double dt);

    //! write results in file
    virtual void WriteResultsInFile();
  

  protected:
  
    Integer size_;        //!< total number of unknowns (equations)


  private:

    // defines subtype of mechanic PDE: plainStrain, 3d, ...
    std::string subType_;

    // Blbablubb
    Integer GetBCDof (const std::string dofStartString);

  };
}

#endif

