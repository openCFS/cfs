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

  //! This class holds information about Coupling terms, such as coupling quantity, values coupling nodes/elements ...
  /*! This class holds information about Coupling terms, such as coupling quantity, values coupling nodes/elements ...
   */
class PDECoupling
{

  // structure for coupling terms
  struct CouplingInterface{
    std::string region;               //!< name of coupling region
    CouplingRegionType regionType ;   //!< type of coupling region (defined in 'environment.hh')
    Integer level;                    //!< multigrid level
    std::vector<Integer> nodes;       //!< vector of coupling nodes 
    std::vector<Elem*> elements;      //!< vector of coupling elements
    Array<Double> values;             //!< array containing coupling values
    Array<Double> oldValues;          //!< array containing couplung values of previous iteration step
    ShortInt dim;                     //!< dimension of coupling values
    Integer size;                     //!< number of coupling values
    NormType normtype;                //!< type of norm
    Double epsilon;                   //!< maximal error tolerance from one step to another
  };

public:
  
  //! constructor
  PDECoupling(Grid * aptgrid, BCs * aptBCs);

  //! destructor
  virtual ~PDECoupling();

  //! register coupling input (used for definition)
  /*!
    \param InType (input) type of input coupling
    \param, Quantity (input) name of input coupling quantity (e.g. 'elecforce')
  */
  virtual void RegisterInput(CouplingInputType InType, 
			     std::string quantity);

  //! add coupling input
  /*!
    \param Quantity (input) name of input coupling quantity
    \param region (input) name of input coupling region
    \param RegionType (input) type of input coupling region
    \param level (input) multigrid level
    \param Couplings (input) vector with all Couplingobjects of coupledpde
  */
  virtual void AddInput(std::string quantity, 
			std::string region,
			CouplingRegionType regionType, 
			Integer level,
			std::vector<PDECoupling*> & couplings);
  
  //! set PDE
  virtual void SetPDE(BasePDE * aPDE);

  //! set coupling output dimension
  virtual void SetOutputDim(Integer i, ShortInt dim);

  //! set coupling output size
  virtual void SetOutputSize(Integer i, Integer size);
    
  //! get PDE name
  virtual std::string GetPDEName();

  //! get number of input couplings
  virtual Integer GetNumInputCouplings();

  //! get number of output couplings
  virtual Integer GetNumOutputCouplings();

  // ------------ input coupling -----------

  //! get input coupling type
  virtual CouplingInputType GetInputType(Integer i)
  { return inputTypes_[i];}

  //! get input coupling quantity
  virtual std::string GetInputQuantity(Integer i)
  { return inputQuantities_[i]; }

  //! get input coupling region
  virtual std::string GetInputRegion(Integer i)
  { return inputInterfaces_[i]->region; }

  //! get input coupling regiontype
  virtual CouplingRegionType GetInputRegionType(Integer i)
  { return inputInterfaces_[i]->regionType; }

  //! get input coupling multigrid level
  virtual Integer GetInputRegionLevel(Integer i)
  { return inputInterfaces_[i]->level; }

 //! get input coupling region nodes
  virtual void GetInputNodes(Integer i, std::vector<Integer>* &Nodes)
  { Nodes  = &(inputInterfaces_[i]->nodes);}

  //! get input coupling region elements
  virtual void GetInputElements(Integer i, std::vector<Elem *>*  &Elements)
  { Elements = &(inputInterfaces_[i]->elements);}

  //! get input coupling values
  virtual void GetInputValues(Integer i, Array<Double>* &Values)
  { Values = &(inputInterfaces_[i]->values);}

  //! get input coupling values
  virtual void GetInputOldValues(Integer i, Array<Double>* &Values)
  { Values = &(inputInterfaces_[i]->oldValues);}

  //! get input coupling values dimension
  virtual ShortInt GetInputDim(Integer i)
  { return inputInterfaces_[i]->dim; }

  //! get input coupling values size
  virtual Integer GetInputSize(Integer i)
  { return inputInterfaces_[i]->size; }

  //! get input coupling norm type
  virtual NormType GetInputNormType(Integer i)
  { return inputInterfaces_[i]->normtype; }

 //! get input coupling epsilon
  virtual Double GetInputEpsilon(Integer i)
  { return inputInterfaces_[i]->epsilon; }

  // ------------- output coupling -----------
  
  //! get output coupling type
  virtual CouplingOutputType GetOutputType(Integer i)
  { return outputTypes_[i];}

  //! get output coupling quantity
  virtual std::string GetOutputQuantity(Integer i)
  { return outputQuantities_[i]; }

  //! get output coupling region
  virtual std::string GetOutputRegion(Integer i)
  { return outputInterfaces_[i]->region; }

  //! get output coupling regiontype
  virtual CouplingRegionType GetOutputRegionType(Integer i)
  { return outputInterfaces_[i]->regionType; }
 
  //! get output coupling multigrid level
  virtual Integer GetOutputRegionLevel(Integer i)
  { return outputInterfaces_[i]->level; }
  
  //! get output coupling region nodes
  virtual void GetOutputNodes(Integer i, std::vector<Integer>* &Nodes)
  { Nodes  = &(outputInterfaces_[i]->nodes);}

  //! get output coupling region elements
  virtual void GetOutputElements(Integer i, std::vector<Elem *>* &Elements)
  { Elements = &(outputInterfaces_[i]->elements);}

  //! get output coupling values
  virtual void GetOutputValues(Integer i, Array<Double>* &Values)
  { Values = &(outputInterfaces_[i]->values);}

  //! get output coupling values
  virtual void GetOutputOldValues(Integer i, Array<Double>* &Values)
  { Values = &(outputInterfaces_[i]->oldValues);}

  //! get output coupling values dimension
  virtual ShortInt GetOutputDim(Integer i)
  { return outputInterfaces_[i]->dim; }

  //! get output coupling values size
  virtual Integer GetOutputSize(Integer i)
  { return outputInterfaces_[i]->size; }
 //! get input coupling norm type
  
  virtual NormType GetOutputNormType(Integer i)
  { return outputInterfaces_[i]->normtype; }

 //! get input coupling epsilon
  virtual Double GetOutputEpsilon(Integer i)
  { return outputInterfaces_[i]->epsilon; }
  
protected:
  
  //! add coupling output
  /*!
    \param OutputType (input) type of output coupling
    \param Quantity (input) name of output coupling quantity
    \param region (input) name of output coupling region
    \param RegionType (input) type of input coupling region
    \param level (input) multigrid level
  */
  virtual CouplingInterface* AddOutput(CouplingOutputType outputType, 
				       std::string quantity, 
				       std::string region,
				       CouplingRegionType regionType,
				       Integer level);

  //! maps input coupling types to output coupling types
  virtual CouplingOutputType Input2OutputType(CouplingInputType inputType);

  BasePDE * myPDE_;                                  //!< pointer to PDE
  Grid * ptGrid_;                                    //!< pointer to grid
  BCs * ptBCs_;                                      //!< pointer to BCs


  // Coupling Output parameters
  std::vector<CouplingOutputType> outputTypes_;      //!< vector containing types of coupling output
  std::vector<std::string> outputQuantities_;        //!< vector containing quantities of coupling output
  std::vector<CouplingInterface*> outputInterfaces_; //!< vector containing pointer to coupling interfaces

  // Coupling Input parameters
  std::vector<CouplingInputType> inputTypes_;        //!< vector containing types of coupling input
  std::vector<std::string> inputQuantities_;         //!< vector conatining quantities of coupling input
  std::vector<CouplingInterface*> inputInterfaces_;  //!< vector containing pointer to coupling interfaces

  // defautl values for coupling
  Double defaultEpsilon;
  NormType defaultNormType;
  

private:

}; // end of class declaration

} // end of namespace

#endif
