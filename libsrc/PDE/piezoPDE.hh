#ifndef FILE_PIEZOPDE
#define FILE_PIEZOPDE


#include "SinglePDE.hh"

#include "Driver/stdSolveStep.hh"

 
namespace CoupledField
{

  //! Class for piezoelectric simulations

  //! Detailed description to follow ...
  class PiezoPDE: public SinglePDE
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

    //! define the SoltionStep-Driver
    virtual void DefineSolveStep();

    /// returns a stiffness integrator appropriate to the actual problem (e.g. 3D)
    BaseForm * GetStiffIntegrator(MaterialData& actSDMat, Boolean reducedInt=FALSE, 
				  Boolean isdamping=FALSE);

    //! do PostProcessing step
    virtual void PostProcess(const Integer level);

    //! write results in file
    //! \param stepOffset offset for starting (time)step
    //! \param timeOffset offset for starting time  
    virtual void WriteResultsInFile(const Integer kstep = 0,
				    const Double asteptime = 0.0,
				    Integer stepOffset = 0,
				    Double timeOffset = 0.0);

    
    // ======================================================
    // COUPLING SECTION
    // ======================================================

    //! initalize PDE coupling
    void InitCoupling(PDECoupling * Coupling)
    { Error ("Coupling not implemented" );}

    //! calculate coupling terms
    void CalcOutputCoupling()
    { Error ("Coupling not implemented" );}

    //! returns if PDE can compute the quantity
    Boolean HasOutput(SolutionType output)
    { Error ("Coupling not implemented" );}
  

   ElemStoreSol<Complex>  GetComplexValuedCharge(){return chargesComplex_;};
    
  protected:
  
    Integer size_;        //!< total number of unknowns (equations)

    //! Obtain information on desired output quantities from parameter file
    
    //! This method is used to query the parameter handling object for the
    //! desired output quantities and translate their literal description into
    //! the internal format by setting the corresponding class attributes.
    //! The output quantities currently supported by the mechanics PDE are
    //! given in the following table. Here 'Keyword' and 'Result Type' refer
    //! to the XML parameter file, while 'Class Attribute' refers to the
    //! internal attribute of the MechPDE class that is set, if the keyword
    //! is specified.\n\n
    //! <table border="1">
    //!   <tr>
    //!     <td><b>Keyword</b></td>
    //!     <td><b>Result Type</b></td>
    //!     <td><b>Class Attribute</b></td>
    //!   </tr>
    //!   <tr>
    //!     <td>displacement</td>
    //!     <td>nodeResults</td>
    //!     <td>savesol_</td>
    //!   </tr>
    //!   <tr>
    //!     <td>velocity</td>
    //!     <td>nodeResults</td>
    //!     <td>savederiv_</td>
    //!   </tr>
    //!   <tr>
    //!     <td>acceleration</td>
    //!     <td>nodeResults</td>
    //!     <td>savederiv2_</td>
    //!   </tr>
    //! </table>
    void ReadStoreResults();

 //! Init the time stepping
    void InitTimeStepping();


  private:

    //postprocessing
    StdVector<std::string> calcEfield_;  //!< contains the subdomains, on which the electric field is computed
    ElemStoreSol<Double> Efield_;  //!< conatins electric field
    ElemStoreSol<Complex> EfieldComplex_;


    StdVector<std::string> calcStress_;  //!< contains the subdomains, on which the stress is computed
    ElemStoreSol<Double> stress_;  //!< conatins mechanical stresses
    ElemStoreSol<Complex> stressComplex_;  //!< conatins mechanical stresses

    ElemStoreSol<Double> charges_;
    ElemStoreSol<Complex> chargesComplex_;
    StdVector<std::string> chargeNeighborRegion_;
    StdVector<std::string> calcCharge_;

    //! contains mechanic velocity
    NodeStoreSol<Double> solDeriv1_;
    
    //! contains mechanic acceleration
    NodeStoreSol<Double> solDeriv2_;
    
    //! calculate Electric field
    void CalcEfield();

    //! calculate complex valued Electric field
    void CalcComplexValuedEfield();
    
    //! calculate stresses
    void CalcStress();

    //! calculate complex valued stresses
    void CalcComplexValuedStress();

    //! calculate Charges
    void CalcCharges();

    //! calculate comlex valued Charges
    void CalcComplexValuedCharges();
    
    
  };
}

#endif

