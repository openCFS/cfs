#ifndef FILE_PDECONF_2003
#define FILE_PDECONF_2003

#include <map>
#include "CoupledPDE/pdecoupling.hh"

#include "PDE/StdPDE.hh"
#include "Utils/StdVector.hh"



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
  virtual void CreateCoupling(StdVector<StdPDE*> &  OrderedPDEs, 
			      StdVector<PDECoupling*> & Couplings,
			      StdVector<StdPDE*> & UnorderedPDEs);

protected:

  //! Defines the possible orderings and coupling possibilities of PDEs
  //! Note: Explicit definition in 'coupledPDE.conf'
  virtual void DefineOrdering();

  Grid * ptGrid_;                        //!< pointer to grid
  BCs * ptBCs_;                          //!< pointer to BC object
  StdVector<Definition*> CoupledPDEs_; //!< vector with coupling definitons

private:
}; // end of class definition
  
  
  
  
//! Helper class (of CoupledPDEDef) for defining ordering of PDEs and coupling terms
class Definition{

public:
  
  //! constructor
  Definition();
  
  //! destructor
  virtual ~Definition();
  
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
				SolutionType Quantity,
				Boolean optional=FALSE);
  
  //! return number of pdes
  virtual Integer GetNumPDEs(){return NumPDEs_;}

  //! return vector with name of PDEs
  /*!
    \param NamePDEs (output) vector containg names of PDEs
  */
  virtual void GetNamePDEs(StdVector<std::string> & NamePDEs)
  {NamePDEs = PDEs_;}
  
  //! return vector with coupling types
   /*!
    \param PDEName (input) name of PDE
    \param InputType (output) vector containing input coupling types of PDE
  */
  virtual void GetCouplingType(std::string  PDEName, 
			       StdVector<CouplingInputType> & InputType) 
  {InputType = InputCouplingTypes_[PDEName];}

  //! return vector with coupling quantities
   /*!
    \param PDEName (input) name of PDE
    \param InputType (output) vector containing input coupling quantities of PDE
  */
  virtual void GetCouplingQuantity(std::string  PDEName, 
				   StdVector<SolutionType> & InputQuantities)
  {InputQuantities = InputCouplingQuantities_[PDEName];}
  

  //! return vector with coupling types
   /*!
    \param PDEName (input) name of PDE
    \param InputType (output) vector containing input coupling types of PDE
  */
  virtual void GetCouplingOptionality(std::string  PDEName, 
				      StdVector<Boolean> & optionality) 
  {optionality = optionalCoupling_[PDEName];}
  
protected:
  //!< number of PDEs
  Integer NumPDEs_;
  
  //!< vector conatining names of PDEs 
  StdVector< std::string> PDEs_;                                             
  
  //!< mapping of PDEnames to input coupling types
  std::map <std::string, StdVector<CouplingInputType> > InputCouplingTypes_;
  
  //!< mapping of PDEnames to input coupling quantities
  std::map <std::string, StdVector<SolutionType> > InputCouplingQuantities_;  
  
  //!< mapping of PDEnames to optionality of input coupling quantities
  std::map <std::string, StdVector<Boolean> > optionalCoupling_;             

};  
} // end of namespace
#endif
