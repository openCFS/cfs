#ifndef FILE_MAGNETICPDE
#define FILE_MAGNETICPDE

#include "General/environment.hh"
#include "basePDE.hh" 

namespace CoupledField
{

  // Forward declaration of classes
  class Coil;

  //! Class for magnetic equation (no adaptivity)
  /*! 
    This class is derived from class BasePDE. 
  */

class MagPDE : public BasePDE
{
public:

  MagPDE(Grid * aptgrid, BCs *aptbcs, TimeFunc *aptTimeFunc,
	 FileType *aptFileType, WriteResults *aptOut);

  //! Deconstructor
  virtual ~MagPDE(){};



  //! define all (bilinearform) integrators needed for this pde
  virtual void DefineIntegrators(const Integer level);

   //! return size of solution
  virtual Integer getSize() const 
  { return numPDENodes_*dofspernode_;}

// ======================================================
// SOLUTION SECTION
// ======================================================

  //! prepare for correct time stepping
  /*! \param dt time step  */
  virtual void InitTimeStepping(const Double dt);

  //!
  virtual void PreStepStatic(const Integer level);

  //! nonlinear static step
  virtual void StepStaticNonLin(const Integer kstep, const Double aTime,
				const Integer level, const Boolean reset);
  //!
  virtual void PostStepStatic(const Integer level);


// ======================================================
// POSTPROCESSING SECTION
// ======================================================

  //!  return pointer to vector with first derivative of solution
  virtual const Vector<Double>& getS1() const { return TS_alg_->GetDeriv1();}

  //! do PostProcessing step
  virtual void PostProcess(const Integer level);

  //! write results in file
  virtual void WriteResultsInFile();

  //! computes the electric energy for each subdomain
  void CalcEnergy();


  //! callculates nodal forces
  void CalcNodeForce(StoreSol<Double> & force, 
		     std::vector<Integer> & nodes, 
		     std::vector<Elem*> & elems,
		     std::vector<std::vector<ShortInt> > & isBoundaryNode,
		     std::vector<std::vector<Integer> > & elemNodeToCouplingNode)
  {Error("CalcNodeForce not implemented",__FILE__,__LINE__);}


  //! GET SOLUTION AT ALL NODES OF AN ELEMENT
  void GetSolOfElement( Vector<Double>& elpot, Vector<Integer>& connect_PDE);

  //! GET 1st derivativ of SOLUTION AT ALL NODES OF AN ELEMENT
  void GetSolDerivOfElement( Vector<Double>& elpot, Vector<Integer>& connect_PDE);

// ======================================================
// COUPLING SECTION
// ======================================================


  //! initalize PDE coupling
  virtual void InitCoupling(PDECoupling * Coupling)
  {Error("InitCoupling not implemented",__FILE__,__LINE__);}
  
  
  //! calculate coupling terms
  virtual void CalcOutputCoupling()
  {Error("CalcOutputCoupling not implemented",__FILE__,__LINE__);}

  //! returns if PDE can compute the quantity
  virtual Boolean HasOutput(std::string output)
  {
    Error("HasOutput not implemented",__FILE__,__LINE__);
    return FALSE;
  }




protected:

  //! Query parameter object for information on coils
  void ReadCoils();

  //! Query parameter object for information on permanent magnets
  void ReadMagnets();

  //!
  void ComputeUI(Vector<Double>& uiSD);
  
  void WriteUI2File(Vector<Double>& uiSD);

  //computes linear part of RHS
  Double SetLinRHS(const Integer level);

  //
  void StoreAlgsysToVec(Vector<Double>& vec, Double * pt);

  /// calculates L2-norm of RHS regarding entries due to penalty formulation
  Double RhsL2Norm(Vector<Double>& stdVec);

  Vector<Double> RhsLinVal_;

  StoreSol<Double> B_;  //!< conatins magnetic field
  StoreSol<Double> Jeddy_;  //!< conatins eddy currents field
  
  // ---- Electric Force variables ---
  StoreSol<Double> Force_;        //!< stores Magnetic force of each element
  std::vector<std::vector<Elem*> > F_Interface_; //!<vector of vectors conaining Elements with acting force
  std::vector<std::vector<std::vector<ShortInt> > > isBoundaryNode_; //!< vector containing flag array for element boundary nodes
  std::vector<std::vector<std::vector<Integer> > > elemNodeToCouplingNode_; //!< assigns each coupling element node the according Coupling Node number
  std::vector<std::vector<Integer> > numBoundaryNodes_;               //!< contains number of surface nodes per element

  // ==========================================================================
  //   COILS
  // ==========================================================================

  //@{ \name Attributes related to coils

  //! Names of coils resp. their subdomains
  std::vector<std::string> coilName_;  

#ifndef XMLPARAMS

  //! vector of parameters describing coils
  std::vector<struct coilDefStruct> coilDef_; 

#else

  //! Parameters of the individual coils;
  std::vector<Coil*> coilDef_;

#endif

  //@}

  // ==========================================================================
  //   PERMANENT MAGNETS
  // ==========================================================================

  //@{ \name Attributes related to permanent magnets

  //! Subdomains containing permanent magnets
  std::vector <std::string> magnetsDomain_;

  //! x-component of direction of magnetisation for each magnet

  //! x-component of direction of magnetisation for each magnet
  //! \todo As suggested by Fred Hofer, the direction of magnetisation of a
  //! permanent magnet must now be specified in the XML parameter file and
  //! no longer in the material data file. While magneticPDE already reads
  //! these data, they are not yet used in the simulation.
  std::vector<Double> magnetsOriX_;

  //! y-component of direction of magnetisation for each magnet
  std::vector<Double> magnetsOriY_;

  //! z-component of direction of magnetisation for each magnet
  std::vector<Double> magnetsOriZ_;

  //@}

  //postprocessing
  std::vector<std::string> calcBfield_;  //!< contains the subdomains, on which the magnetic field is computed
  std::vector<std::string> calcEnergy_;  //!< contains the subdomains, on which the magnetic energy is computed
  std::vector<std::string> calcEddy_;  //!< contains the subdomains, on which the eddy currents are computed

  std::ofstream * UIfile_; //!< file for informational output
  std::string UIfilename_;      //!< name of file for saving current/voltage values
   
  private:

#ifdef XMLPARAMS
    //! Obtain information on desired output quantities from parameter file

    //! This method is used to query the parameter handling object for the
    //! desired output quantities and translate their literal description into
    //! the internal format by setting the corresponding class attributes.
    //! The output quantities currently supported by the electrostatics PDE are
    //! given in the following table. Here 'Keyword' and 'Result Type' refer
    //! to the XML parameter file, while 'Class Attribute' refers to the
    //! internal attribute of the MagPDE class that is set, if the keyword
    //! is specified.\n\n
    //! <table border="1">
    //!   <tr>
    //!     <td><b>Physical quantity</b></td>
    //!     <td><b>Keyword in parameter file</b></td>
    //!     <td><b>Result Type</b></td>
    //!     <td><b>Class Attribute</b></td>
    //!   </tr>
    //!   <tr>
    //!     <td>Magnetic flux density \f$\vec{B}\f$</td>
    //!     <td>bfield</td>
    //!     <td>Element result</td>
    //!     <td>calcBfield_</td>
    //!   </tr>
    //!   <tr>
    //!     <td>Magnetic energy \f$W\f$</td>
    //!     <td>energy</td>
    //!     <td>Element result</td>
    //!     <td>calcEnergy_</td>
    //!   </tr>
    //!   <tr>
    //!     <td>Eddy current density \f$I_E\f$</td>
    //!     <td>eddy</td>
    //!     <td>Element result</td>
    //!     <td>calcEddy_</td>
    //!   </tr>
    //!   <tr>
    //!     <td>depends on formulation</td>
    //!     <td>none -- always true</td>
    //!     <td>Nodal results</td>
    //!     <td>savesol_</td>
    //!   </tr>
    //!   <tr>
    //!     <td align="center">---</td>
    //!     <td>none -- always false</td>
    //!     <td>Nodal results</td>
    //!     <td>deriv1_</td>
    //!   </tr>
    //!   <tr>
    //!     <td align="center">---</td>
    //!     <td>none -- always false</td>
    //!     <td>Nodal results</td>
    //!     <td>deriv2_</td>
    //!   </tr>
    //! </table>
    void ReadStoreResults();
#endif

};

} // end of namespace
#endif

