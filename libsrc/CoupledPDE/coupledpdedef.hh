#ifndef FILE_PDECONF_2003
#define FILE_PDECONF_2003

#include <vector>
#include <map>
#include <CoupledPDE/pdecoupling.hh>

#ifndef NEWBASEPDE
#include <PDE/basepde.hh>
#else
#include <PDE/newBasePDE.hh>
#endif //#ifndef NEWBASEPDE



namespace CoupledField
{

class Definition;
class Grid;
class BCs;


//! Class for composing coupled PDEs out of single field PDEs
/*! This class arranges unordered PDEs in the right order and defines the coupling possibilities
 */
class CoupledPDEDef
{
public:
  
  //! Constructor
  CoupledPDEDef(Grid * aptGrid, BCs * aptBCs);

  //! Destructor
  virtual ~CoupledPDEDef();

  //! Creates coupling between PDEs
  /*!
    \param OrderedPDEs (output) vector of PDEs in right ordering for solving
    \param Couplings (output) vector with coupling objects
    \param UnorderedPDEs (inpu) unordered vector with PDEs for coupling
  */
  virtual void CreateCoupling(std::vector<BasePDE*> &  OrderedPDEs, 
			      std::vector<PDECoupling*> & Couplings,
			      std::vector<BasePDE*> & UnorderedPDEs);

protected:

  //! Defines the possible orderings and coupling possibilities of PDEs
  //! Note: Explicit definition in 'coupledPDE.conf'
  virtual void DefineOrdering();

  Grid * ptGrid_;                        //!< pointer to grid
  BCs * ptBCs_;                          //!< pointer to BC object
  std::vector<Definition*> CoupledPDEs_; //!< vector with coupling definitons

private:
}; // end of class definition
  
  
  
  
//! Helper class (of CoupledPDEDef) for defining ordering of PDEs and coupling terms
class Definition{

public:
  
  //! constructor
  Definition();
  
  //! destructor
  ~Definition();
  
  //! adds a PDE to definition
  virtual void AddPDE(std::string);
  
  //! adds input coupling quantities
  /*!
    \param PDEName (input) name of PDE
    \param InType (input) type of input coupluing (coordinates, RHS, ..)
    \param Quantity (input) name of input coupling quantity (e.g. 'elecforce')
  */
  virtual void AddInputCoupling(std::string PDEName, 
				CouplingInputType InType, 
				std::string Quantity);
  
  //! return number of pdes
  virtual Integer GetNumPDEs(){return NumPDEs_;}

  //! return vector with name of PDEs
  /*!
    \param NamePDEs (output) vector containg names of PDEs
  */
  virtual void GetNamePDEs(std::vector<std::string> & NamePDEs)
  {NamePDEs = PDEs_;}
  
  //! return vector with coupling types
   /*!
    \param PDEName (input) name of PDE
    \param InputType (output) vector containing input coupling types of PDE
  */
  virtual void GetCouplingType(std::string  PDEName, 
			       std::vector<CouplingInputType> & InputType) 
  {InputType = InputCouplingTypes_[PDEName];}

  //! return vector with coupling quantities
   /*!
    \param PDEName (input) name of PDE
    \param InputType (output) vector containing input coupling quantities of PDE
  */
  virtual void GetCouplingQuantity(std::string  PDEName, 
				   std::vector<std::string> & InputQuantities)
  {InputQuantities = InputCouplingQuantities_[PDEName];}
  
  
protected:
  Integer NumPDEs_;                                                            //!< number of PDEs
  std::vector< std::string> PDEs_;                                             //!< vector conatining names of PDEs 
  std::map <std::string, std::vector<CouplingInputType> > InputCouplingTypes_; //!< mapping of PDEnames to input coupling types
  std::map <std::string, std::vector<std::string> > InputCouplingQuantities_;  //!< mapping of PDEnames to input coupling quantities

};  
} // end of namespace
#endif
