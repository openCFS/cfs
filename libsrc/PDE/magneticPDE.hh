#ifdef NEWBASEPDE


#ifndef FILE_MAGNETICPDE
#define FILE_MAGNETICPDE

#include "General/environment.hh"
#include "newBasePDE.hh" 

namespace CoupledField
{

  //! Class for magnetic equation (no adaptivity)
  /*! 
    This class is derived from class BasePDE. 
  */

class MagPDE : public BasePDE
{
public:

  //! Constructor. here we read integration parameters
  /*!
    \param 
    \param aGrid pointer to grid
    \param aBCs pointer to Boundary condition object
    \param aGrid pointer to class Grid
    \param aInFile pointer to class FileType. input data.Boolean MagPDE::HasOutput(std::string output)
{
#ifdef TRACE
  (*trace) << "entering MagPDE::HasOutput" << std::endl;
#endif
  
  if (output == "elecforce")
    return TRUE;

  if (output == "elecpotential")
    return TRUE;

  if (output == "elecfield")
    return TRUE;

  return FALSE;

    \param aOutFile  pointer to class WriteResults. output data.
    \param aTimeFunc pointer to class TimeFunc
  */
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

  //!
  virtual void PreStepStatic(const Integer level);

  //!
  virtual void PostStepStatic(const Integer level);


// ======================================================
// POSTPROCESSING SECTION
// ======================================================


  //! do PostProcessing step
  virtual void PostProcess(const Integer level);

  //! write results in file
  virtual void WriteResultsInFile();

  //! computes the electric energy for each subdomain
  void CalcEnergy();


  //! callculates nodal forces
  void CalcNodeForce(Array<Double> & force, 
		     std::vector<Integer> & nodes, 
		     std::vector<Elem*> & elems,
		     std::vector<std::vector<ShortInt> > & isBoundaryNode,
		     std::vector<std::vector<Integer> > & elemNodeToCouplingNode)
  {Error("CalcNodeForce not implemented");}


  //! GET SOLUTION AT ALL NODES OF AN ELEMENT
  void GetSolOfElement( Vector<Double>& elpot, Vector<Integer>& connect_PDE);


// ======================================================
// COUPLING SECTION
// ======================================================


  //! initalize PDE coupling
  virtual void InitCoupling(PDECoupling * Coupling)
  {Error("InitCoupling not implemented");}
  
  
  //! calculate coupling terms
  virtual void CalcOutputCoupling()
  {Error("CalcOutputCoupling not implemented");}

  //! returns if PDE can compute the quantity
  virtual Boolean HasOutput(std::string output)
  {Error("HasOutput not implemented");}




protected:

  //! reads all data in the config-file belonging to coils
  void ReadCoils();

  Array<Double> B_;  //!< conatins magnetic field
  
  // ---- Electric Force variables ---
  Array<Double> Force_;        //!< stores Magnetic force of each element
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
};

} // end of namespace
#endif

#endif //#ifdef NEWBASEPDE
