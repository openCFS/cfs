#ifndef FILE_MAGNETICPDE
#define FILE_MAGNETICPDE

#include "General/environment.hh"
#include "basePDE.hh" 

namespace CoupledField
{

  //! Class for magnetic equation (no adaptivity)
  /*! 
    This class is derived from class BasePDE. 
  */

class MagPDE : public BasePDE
{
public:

  MagPDE(Grid * aptgrid, BCs *aptbcs, TimeFunc *aptTimeFunc, FileType *aptFileType, WriteResults *aptOut);

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

  //! reads all data in the config-file belonging to coils
  void ReadCoils();

  //!
  void ComputeUI(Vector<Double>& uiSD);
  
  void WriteUI2File(Vector<Double>& uiSD);


  StoreSol<Double> B_;  //!< conatins magnetic field
  StoreSol<Double> Jeddy_;  //!< conatins eddy currents field
  
  // ---- Electric Force variables ---
  StoreSol<Double> Force_;        //!< stores Magnetic force of each element
  std::vector<std::vector<Elem*> > F_Interface_; //!<vector of vectors conaining Elements with acting force
  std::vector<std::vector<std::vector<ShortInt> > > isBoundaryNode_; //!< vector containing flag array for element boundary nodes
  std::vector<std::vector<std::vector<Integer> > > elemNodeToCouplingNode_; //!< assigns each coupling element node the according Coupling Node number
  std::vector<std::vector<Integer> > numBoundaryNodes_;               //!< contains number of surface nodes per element

  // coils
  std::vector <std::string> coilDomain_;  //!< name of all subdomains containing coils
  std::vector<struct coilDefStruct> coilDef_; //!< vector of paramters describing coils

  // permanent magnets
  std::vector <std::string> magnetsDomain_;  //!< name of all subdomains containing permanent magnets

  //postprocessing
  std::vector<std::string> calcBfield_;  //!< contains the subdomains, on which the magnetic field is computed
  std::vector<std::string> calcEnergy_;  //!< contains the subdomains, on which the magnetic energy is computed
  std::vector<std::string> calcEddy_;  //!< contains the subdomains, on which the eddy currents are computed

  std::ofstream * UIfile_; //!< file for informational output
  std::string UIfilename_;      //!< name of file for saving current/voltage values
   
};

} // end of namespace
#endif

