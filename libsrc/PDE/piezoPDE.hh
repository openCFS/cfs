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
    BaseForm * GetStiffIntegrator(MaterialData& actSDMat, Boolean reducedInt=FALSE, 
				  Boolean isdamping=FALSE);

    //! Init the time stepping
    //! \param dt time step
    virtual void InitTimeStepping(const Double dt);
    
    //! return pointer to vector with second derivative of solution
    virtual const Vector<Double> & getS2() const { return TS_alg_->GetDeriv2();}

    //! do PostProcessing step
    virtual void PostProcess(const Integer level);

    //! write results in file
    //! \param stepOffset offset for starting (time)step
    //! \param timeOffset offset for starting time  
    virtual void WriteResultsInFile(Integer stepOffset = 0,
				    Double timeOffset = 0.0);
    
  protected:
  
    Integer size_;        //!< total number of unknowns (equations)

#ifdef XMLPARAMS
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
#endif
  private:

    //postprocessing
    StdVector<std::string> calcEfield_;  //!< contains the subdomains, on which the electric field is computed
    ElemStoreSol<Double> Efield_;  //!< conatins electric field


    StdVector<std::string> calcStress_;  //!< contains the subdomains, on which the stress is computed
    ElemStoreSol<Double> stress_;  //!< conatins mechanical stresses

    ElemStoreSol<Double> charges_;
    StdVector<std::string> chargeNeighborRegion_;
    StdVector<std::string> calcCharge_;

    //! contains mechanic velocity
    NodeStoreSol<Double> solDeriv1_;
    
    //! contains mechanic acceleration
    NodeStoreSol<Double> solDeriv2_;
    
    //! calculate stresses
    void CalcEfield();
    
    //! calculate stresses
    void CalcStress();

    //! calculate stresses
    void CalcCharges();
    
    StdVector<std::string> pressSurf_;  //!< surface of pressure loads
    StdVector<Double>      pressVals_;  //!< values of the pressure loads
    StdVector<std::string> pressFnc_;   //!< function names of pressure loads
    
  };
}

#endif

