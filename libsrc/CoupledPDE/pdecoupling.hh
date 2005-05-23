#ifndef FILE_PDECOUPLING_2003
#define FILE_PDECOUPLING_2003

#include <list>

#include "General/environment.hh"
#include "Utils/vector.hh"
#include "Utils/StdVector.hh"
#include "DataInOut/MaterialData.hh"

namespace CoupledField
{

  // Forward declarations of classes
  class StdPDE;
  class Grid;
  class Elem;
  class CouplingMemento;
  template<class TYPE> class Vector;
  template<class TYPE> class Matrix;

  //! This class holds information about Coupling terms, such as coupling quantity, 
  //! values coupling nodes/elements ...

  //! This class holds information about Coupling terms, such as coupling quantity, 
  //! values coupling nodes/elements ...
  class PDECoupling
  {

    // friend declaration
    friend class CouplingMemento;

    // structure for coupling terms
    struct CouplingInterface{
    public:
    
      //! constructor
      CouplingInterface();
    
      //! copy constructor
      CouplingInterface(const CouplingInterface & x);

      //! destructor
      ~CouplingInterface();

      //! assignement operator
      CouplingInterface & operator=(const CouplingInterface & x);

      //! name of coupling region
      StdVector<std::string> regions;       

      //! type of coupling region (defined in 'environment.hh')
      CouplingRegionType regionType;        

      //! vector of coupling nodes 
      StdVector<Integer> nodes;           

      //! vector of coupling elements
      StdVector<Elem*> elements;          

      //! vector containing neighbouring region of interface
      //! w.r.t. pde which calculates the values for the coupling
      StdVector<std::string> neighInputRegions; 
    
      //! vector of materials at coupling interface
      StdVector<MaterialData*> materials; 

      //! vector with materials for opposite PDE
      //! ( = pde, which uses this interface as input )
      StdVector<MaterialData*>  oppositePdeMaterials;         

      //! array containing coupling values
      CFSVector * values;

      //! array containing coupling values of previous iteration step
      CFSVector * oldValues;

      //! dof of coupling values
      ShortInt dof;

      //! number of couplingnodes
      Integer numNodes;

      //! number of couplingelements
      Integer numElems;

      //! type of norm
      NormType normtype;

      //! maximal error tolerance from one step to another
      Double epsilon;
    };

  public:
  
    //! constructor
    PDECoupling( Grid * aptgrid );

    //! destructor
    virtual ~PDECoupling();

    //! register coupling input (used for definition)
    /*!
      \param InType (input) type of input coupling
      \param, Quantity (input) name of input coupling quantity (e.g. 'elecforce')
    */
    virtual void RegisterInput(CouplingInputType InType, 
                               SolutionType quantity);

    //! add coupling input
    /*!
      \param numCoupling (input) number of Couplinginterface
      \param Quantity (input) name of input coupling quantity
      \param interfaceNames (input) names of input coupling interfaces
      \param RegionType (input) type of input coupling interfaces
      \param epsilon (input) tolerance for quantity
      \param normtype (input) normtype of epsilon
      \param Couplings (input) vector with all Couplingobjects of coupledpde
    */
    virtual void AddInput(SolutionType quantity, 
                          StdVector<std::string> &interfaceNames,
                          CouplingRegionType regionType, 
                          StdVector<std::string> & neighRegions,
                          Double epsilon,
                          NormType normtype,
                          StdVector<PDECoupling*> & couplings);
  
    //! set PDE
    virtual void SetPDE(StdPDE * aPDE);

    //! set coupling output dof
    virtual void SetOutputDof(Integer i, ShortInt dof);

    //! set output coupling nodes size
    virtual void SetOutputNumNodes(Integer i, Integer size);
 
    //! set output coupling elements size
    virtual void SetOutputNumElems(Integer i, Integer size);  

    //! get PDE name
    virtual std::string GetPDEName();
  
    //! get number of input couplings
    virtual Integer GetNumInputCouplings();

    //! get number of output couplings
    virtual Integer GetNumOutputCouplings();


    //! creates a new CFSVector-object for the 
    //! couplingvalues.
    //! This method has to be called from the according
    //! PDE in the method 'InitCoupling'
    /*!
      \param i (input) Number of Couplinginterface
      \param solType (input) Type of Solution (ref. Enum-type)
      \param isComplex (inut) True if values are complex
    */
    virtual void CreateCouplingVector(Integer i,
                                      Boolean isComplex);
  
    // ------------ input coupling -----------

    //! get input coupling type
    virtual CouplingInputType GetInputType(Integer i)
    { return inputTypes_[i];}

    //! get input coupling quantity
    virtual SolutionType GetInputQuantity(Integer i)
    { return inputQuantities_[i]; }

    //! get input coupling region
    virtual void GetInputRegions(Integer i, StdVector<std::string> & regions)
    { regions = inputInterfaces_[i]->regions; }

    //! get input coupling regiontype
    virtual CouplingRegionType GetInputRegionType(Integer i)
    { return inputInterfaces_[i]->regionType; }

    //! get input coupling region nodes
    virtual void GetInputNodes(Integer i, StdVector<Integer>* &nodes)
    { nodes  = &(inputInterfaces_[i]->nodes);}

    //! get input coupling region elements
    virtual void GetInputElements(Integer i, StdVector<Elem *>*  &elements)
    { elements = &(inputInterfaces_[i]->elements);}

    //! get input neighbour region
    virtual void GetInputNeighbourRegion(Integer i, StdVector<std::string>* &regions)
    { regions = &(inputInterfaces_[i]->neighInputRegions);}

    //! get input coupling values
    virtual void GetInputValues(Integer i, CFSVector* &values)
    { values = (inputInterfaces_[i]->values);}

    //! get input coupling values
    virtual void GetInputOldValues(Integer i, CFSVector* &values)
    { values = (inputInterfaces_[i]->oldValues);}

    //! get input coupling values dof
    virtual ShortInt GetInputDof(Integer i)
    { return inputInterfaces_[i]->dof; }

    //! get input coupling nodes size
    virtual Integer GetInputNumNodes(Integer i)
    { return inputInterfaces_[i]->numNodes; }

    //! get input coupling elems size
    virtual Integer GetInputNumElems(Integer i)
    { return inputInterfaces_[i]->numElems; }

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
    virtual SolutionType GetOutputQuantity(Integer i)
    { return outputQuantities_[i]; }

    //! get output coupling region
    virtual void GetOutputRegions(Integer i, StdVector<std::string> &regions)
    { regions = outputInterfaces_[i]->regions; }

    //! get output coupling regiontype
    virtual CouplingRegionType GetOutputRegionType(Integer i)
    { return outputInterfaces_[i]->regionType; }
 
    //! get output coupling region nodes
    virtual void GetOutputNodes(Integer i, StdVector<Integer>* &nodes)
    { nodes  = &(outputInterfaces_[i]->nodes);}

    //! get output coupling region elements
    virtual void GetOutputElements(Integer i, StdVector<Elem *>* &elements)
    { elements = &(outputInterfaces_[i]->elements);}
  
    //! get output neighbour elements
    //virtual void GetOutputNeighbourElems(Integer i, StdVector<Elem *>*  &elements)
    //{ elements = &(outputInterfaces_[i]->neighbours);}
  
    //! get output neighbour region
    //virtual void GetOutputNeighbourRegion(Integer i, StdVector<std::string>* &regions)
    //{ regions = &(outputInterfaces_[i]->neighInputRegions);}
    
    //! get output coupling region materials
    virtual void GetMaterials(Integer i, StdVector<MaterialData *>* &mat)
    { mat = &(outputInterfaces_[i]->materials);} 
  
    //! get output coupling region materials of opposite pde
    virtual void GetOppositeMaterials(Integer i, StdVector<MaterialData *>* &mat)
    { mat = &(outputInterfaces_[i]->oppositePdeMaterials);}

    //! get output coupling values
    virtual void GetOutputValues(Integer i, CFSVector* &values)
    { values = (outputInterfaces_[i]->values);}

    //! get old output coupling values
    virtual void GetOutputOldValues(Integer i, CFSVector* &values)
    { values = (outputInterfaces_[i]->oldValues);}

    //! get output coupling values dof
    virtual ShortInt GetOutputDof(Integer i)
    { return outputInterfaces_[i]->dof; }

    //! get output coupling number nodes
    virtual Integer GetOutputNumNodes(Integer i)
    { return outputInterfaces_[i]->numNodes; }

    //! get output coupling number elems
    virtual Integer GetOutputNumElems(Integer i)
    { return outputInterfaces_[i]->numElems; }

    //! get output coupling norm type
    virtual NormType GetOutputNormType(Integer i)
    { return outputInterfaces_[i]->normtype; }

    //! get output coupling epsilon
    virtual Double GetOutputEpsilon(Integer i)
    { return outputInterfaces_[i]->epsilon; }

    // ------------- memento operations  -----------

    //! get the encapsulated state of the Coupling object
  
    //! returns the current state of the Coupling object.
    //! This is needed to
    //! enable full MultiSequence simulation, where from one step to 
    //! another the coupling values have to be passed.
    //! The CouplingMemento object encapsulates this information. 
    //! Later on the information can be given back to the Coupling object
    //! with the method SetMemento();
    //! \param memento (output) Object where the current state gets saved
    virtual void GetMemento(CouplingMemento & memento);
  
    //! set the encapsulated state of the Coupling object
  
    //! Set the saved state of the Coupling obejct, which was previously
    //! stored in a CouplingMemento object. This is needed to
    //! enable full MultiSequence simulation, where from one step to 
    //! another the coupling values have to be passed.
    //! The CouplingMemento object encapsulates this information. 
    //! \param memento (input) Previously saved state of the coupling object
    virtual void SetMemento(CouplingMemento & memento); 

  protected:
  
    //! add coupling output
    /*!
      \param OutputType (input) type of output coupling
      \param Quantity (input) name of output coupling quantity
      \param region (input) name of output coupling region
      \param RegionType (input) type of input coupling region
    */
    virtual CouplingInterface* AddOutput(CouplingOutputType outputType, 
                                         SolutionType quantity, 
                                         StdVector<std::string> & region,
                                         CouplingRegionType regionType);
                                      
    //! maps input coupling types to output coupling types
    virtual CouplingOutputType Input2OutputType(CouplingInputType inputType);

    StdPDE * myPDE_;                                  //!< pointer to PDE
    Grid * ptGrid_;                                    //!< pointer to grid


    // Coupling Output parameters
    StdVector<CouplingOutputType> outputTypes_;      //!< vector containing types of coupling output
    StdVector<SolutionType> outputQuantities_;        //!< vector containing quantities of coupling output
    StdVector<CouplingInterface*> outputInterfaces_; //!< vector containing pointer to coupling interfaces

    // Coupling Input parameters
    StdVector<CouplingInputType> inputTypes_;        //!< vector containing types of coupling input
    StdVector<SolutionType> inputQuantities_;         //!< vector conatining quantities of coupling input
    StdVector<CouplingInterface*> inputInterfaces_;  //!< vector containing pointer to coupling interfaces
  
    // defautl values for coupling
    Double defaultEpsilon;
    NormType defaultNormType;
  

  private:

  }; // end of class declaration


} // end of namespace

#endif
