#ifndef FILE_PDECONF_2003
#define FILE_PDECONF_2003

#include <vector>
#include <map>
#include <PDE/basepde.hh>
#include <CoupledPDE/pdecoupling.hh>



namespace CoupledField
{

class Definition;
class Grid;
class BCs;


//! Class for composing coupled PDEs out of sinfle field PDEs
/*! This class arranges unordered PDEs in the right order and defines the coupling possibilities
 */
class CoupledPDEDef
{
public:
  
  //! Constructor
  CoupledPDEDef(Grid * aptGrid, BCs * aptBCs);

  //! Destructor
  virtual ~CoupledPDEDef();

  //! 
  virtual void CreateCoupling(std::vector<BasePDE*> &  OrderedPDEs, 
			      std::vector<PDECoupling*> & Couplings,
			      std::vector<BasePDE*> & UnorderedPDEs);

protected:

  virtual void DefineOrdering();

  Grid * ptGrid_;
  BCs * ptBCs_;

  std::vector<Definition*> CoupledPDEs_;

private:
}; // end of class definition




//! Helper class (of CoupledPDEDef) for defining order of PDEs and coupling terms
class Definition{

public:
  
  //! constructor
  Definition();
  
  //! destructor
  ~Definition();
  
  //! adds a PDE t
  virtual void AddPDE(std::string);

  //! adds input coupling quantities
  virtual void AddInputCoupling(std::string PDEName, CouplingInputType InType, std::string Quantity);

  //! return number of pdes
  virtual Integer GetNumPDEs(){return NumPDEs_;}

  //! return vector with name of PDEs
  virtual void GetNamePDEs(std::vector<std::string> & NamePDEs)
  {NamePDEs = PDEs_;}
  
  //! return vector with coupling types
  virtual void GetCouplingType(std::string  PDEName, std::vector<CouplingInputType> & InputType) 
  {InputType = InputCouplingTypes_[PDEName];}

  //! return vector with coupling quantities
  virtual void GetCouplingQuantity(std::string  PDEName, std::vector<std::string> & InputQuantities)
  {InputQuantities = InputCouplingQuantities_[PDEName];}
  
  
protected:
  Integer NumPDEs_;
  std::vector< std::string> PDEs_;
  std::map <std::string, std::vector<CouplingInputType> > InputCouplingTypes_;
  std::map <std::string, std::vector<std::string> > InputCouplingQuantities_;

};  
} // end of namespace
#endif
