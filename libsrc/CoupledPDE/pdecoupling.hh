#ifndef FILE_PDECOUPLING_2003
#define FILE_PDECOUPLING_2003

#include <General/environment.hh>
#include <Utils/vector.hh>
#include <Utils/array.hh>
#include <list>

namespace CoupledField
{

class BasePDE;
class Grid;
class Elem;
class BCs;


  //! This class holds information about Coupling terms, such as coupling quantity, values
  //! coupling nodes/elements ...
class PDECoupling
{

  // structure for coupling terms
  struct CouplingInterface{
    std::string Region;
    CouplingRegionType RegionType ;
    Integer level;
    std::vector<Integer> Nodes;
    std::vector<Elem*> Elements;
    Array<Double> Values;
    ShortInt Dim;
    Integer Size;
  };

public:
  
  //! constructor
  PDECoupling(Grid * aptgrid, BCs * aptBCs);

  //! destructor
  virtual ~PDECoupling();

  //! register coupling input
  virtual void RegisterInput(CouplingInputType InType, 
			     std::string Quantity);

  //! add coupling input
  virtual void AddInput(std::string Quantity, 
			std::string region,
			CouplingRegionType RegionType, 
			Integer level,
			std::vector<PDECoupling*> & Couplings);
  
  //! set PDE
  virtual void SetPDE(BasePDE * aPDE);

  //! set coupling output dimension
  virtual void SetOutputDim(Integer i, ShortInt Dim);

  //! set coupling output size
  virtual void SetOutputSize(Integer i, Integer Size);
    
  //! get number of input couplings
  virtual Integer GetNumInputCouplings();

  //! get number of output couplings
  virtual Integer GetNumOutputCouplings();

  // ------------ input coupling -----------

  //! get input coupling type
  virtual CouplingInputType GetInputType(Integer i)
  { return InputTypes_[i];}

  //! get input coupling quantity
  virtual std::string GetInputQuantity(Integer i)
  { return InputQuantities_[i]; }

  //! get input coupling region
  virtual std::string GetInputRegion(Integer i)
  { return InputInterfaces_[i]->Region; }

  //! get input coupling regiontype
  virtual CouplingRegionType GetInputRegionType(Integer i)
  { return InputInterfaces_[i]->RegionType; }

  //! get input coupling multigrid level
  virtual Integer GetInputRegionLevel(Integer i)
  { return InputInterfaces_[i]->level; }

 //! get input coupling region nodes
  virtual void GetInputNodes(Integer i, std::vector<Integer>* &Nodes)
  { Nodes  = &(InputInterfaces_[i]->Nodes);}

  //! get input coupling region elements
  virtual void GetInputElements(Integer i, std::vector<Elem *>*  &Elements)
  { Elements = &(InputInterfaces_[i]->Elements);}

  //! get input coupling values
  virtual void GetInputValues(Integer i, Array<Double>* &Values)
  { Values = &(InputInterfaces_[i]->Values);}

  //! get input coupling values dimension
  virtual ShortInt GetInputDim(Integer i)
  { return InputInterfaces_[i]->Dim; }

  //! get input coupling values size
  virtual Integer GetInputSize(Integer i)
  { return InputInterfaces_[i]->Size; }
 
  // ------------- output coupling -----------
  
  //! get output coupling type
  virtual CouplingOutputType GetOutputType(Integer i)
  { return OutputTypes_[i];}

  //! get output coupling quantity
  virtual std::string GetOutputQuantity(Integer i)
  { return OutputQuantities_[i]; }

  //! get output coupling region
  virtual std::string GetOutputRegion(Integer i)
  { return OutputInterfaces_[i]->Region; }

  //! get output coupling regiontype
  virtual CouplingRegionType GetOutputRegionType(Integer i)
  { return OutputInterfaces_[i]->RegionType; }
 
  //! get output coupling multigrid level
  virtual Integer GetOutputRegionLevel(Integer i)
  { return OutputInterfaces_[i]->level; }
  
  //! get output coupling region nodes
  virtual void GetOutputNodes(Integer i, std::vector<Integer>* &Nodes)
  {
    //std::cerr << "In PDECoupling(PDE: " ;
    std::cerr << " Accessing Element " << i;
    std::cerr << " Size of Nodes: " << OutputInterfaces_[i]->Nodes.size() << std::endl;
    Nodes  = &(OutputInterfaces_[i]->Nodes);}

  //! get output coupling region elements
  virtual void GetOutputElements(Integer i, std::vector<Elem *>* &Elements)
  { Elements = &(OutputInterfaces_[i]->Elements);}

  //! get output coupling values
  virtual void GetOutputValues(Integer i, Array<Double>* &Values)
  { Values = &(OutputInterfaces_[i]->Values);}

  //! get output coupling values dimension
  virtual ShortInt GetOutputDim(Integer i)
  { return OutputInterfaces_[i]->Dim; }

  //! get output coupling values size
  virtual Integer GetOutputSize(Integer i)
  { return OutputInterfaces_[i]->Size; }

   //! calculates norm of 
  virtual Double CalcNorm(std::string type);


protected:
  
  //! add coupling output
  virtual CouplingInterface* AddOutput(CouplingOutputType OutputType, 
				       std::string Quantity, 
				       std::string region,
				       CouplingRegionType RegionType,
				       Integer level);

  virtual CouplingOutputType Input2OutputType(CouplingInputType InputType);

  BasePDE * myPDE_;
  Grid * ptGrid_;
  BCs * ptBCs_;


  // Coupling Output parameters
  std::vector<CouplingOutputType> OutputTypes_;     //!< vector containing types of coupling output
  std::vector<std::string> OutputQuantities_;         //!< vector containing quantities of coupling output
  std::vector<CouplingInterface*> OutputInterfaces_; //!< vector containing pointer to coupling interfaces

  // Coupling Input parameters
  std::vector<CouplingInputType> InputTypes_;       //!< vector containing types of coupling input
  std::vector<std::string> InputQuantities_;          //!< vector conatining quantities of coupling input
  std::vector<CouplingInterface*> InputInterfaces_;  //!< vector containing pointer to coupling interfaces


private:

}; // end of class declaration

} // end of namespace

#endif
