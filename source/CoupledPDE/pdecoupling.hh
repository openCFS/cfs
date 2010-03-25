// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_PDECOUPLING_2003
#define FILE_PDECOUPLING_2003

#include "General/environment.hh"
#include "MatVec/vector.hh"
#include "Utils/StdVector.hh"
#include "Materials/baseMaterial.hh"
#include "Utils/boost-serialization.hh"

namespace CoupledField
{

  // Forward declarations of classes
  class StdPDE;
  class Grid;
  class Elem;
  class CouplingMemento;
  template<class TYPE> class Vector;
  template<class TYPE> class Matrix;

   

  //! This class holds information about itertive Coupling terms
  class PDECoupling
  {

    // friend declaration
    friend class CouplingMemento;

  public:

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
      StdVector<UInt> nodes;           

      //! vector of coupling elements
      StdVector<Elem*> elements;          

      //! vector containing neighbouring region of interface
      //! w.r.t. pde which calculates the values for the coupling
      StdVector<std::string> neighInputRegions; 
    
      //! vector of materials at coupling interface
      StdVector<BaseMaterial*> materials; 

      //! vector with materials for opposite PDE
      //! ( = pde, which uses this interface as input )
      StdVector<BaseMaterial*>  oppositePdeMaterials;         

      //! array containing coupling values
      SingleVector * values;

      //! array containing coupling values of previous iteration step
      SingleVector * oldValues;

      //! array containing coupling values of n-1 (previus) time step
      SingleVector * values_tn_1;

      //! array containing coupling values of n-2 (previus) time step
      SingleVector * values_tn_2;

      //! array containing coupling values of n-3 (previus) time step
      SingleVector * values_tn_3;

      //! dof of coupling values
      UInt dof;

      //! number of couplingnodes
      UInt numNodes;

      //! number of couplingelements
      UInt numElems;

      //! type of norm
      NormType normtype;

      //! maximal error tolerance from one step to another
      Double epsilon;
    private:

      // =======================================================================
      // SERIALIZATION FUNCTIONS
      // =======================================================================
      // These functions allow us to write a memento directly
      // into an boost::archive, for saving on a disk or in a 
      // iostream object
      
      //! allow serialization class to access memento entries
      friend class boost::serialization::access;
      
      //! Saving internal state into a boost::archive
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version);
      
    };

  
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
    virtual void SetOutputDof(UInt i, UInt dof);

    //! set output coupling nodes size
    virtual void SetOutputNumNodes(UInt i, UInt size);
 
    //! set output coupling elements size
    virtual void SetOutputNumElems(UInt i, UInt size);  

    //! get pointer to PDE
    virtual StdPDE* GetPDE();
  
    //! get number of input couplings
    virtual UInt GetNumInputCouplings();

    //! get number of output couplings
    virtual UInt GetNumOutputCouplings();


    //! creates a new SingleVector-object for the 
    //! couplingvalues.
    //! This method has to be called from the according
    //! PDE in the method 'InitCoupling'
    /*!
      \param i (input) Number of Couplinginterface
      \param solType (input) Type of Solution (ref. Enum-type)
      \param isComplex (inut) True if values are complex
    */
    virtual void CreateCouplingVector(UInt i,
                                      bool isComplex);
  
    // ------------ input coupling -----------

    //! get input coupling type
    virtual CouplingInputType GetInputType(UInt i)
    { return inputTypes_[i];}

    //! get input coupling quantity
    virtual SolutionType GetInputQuantity(UInt i)
    { return inputQuantities_[i]; }

    //! get input coupling region
    virtual void GetInputRegions(UInt i, StdVector<std::string> & regions)
    { regions = inputInterfaces_[i]->regions; }

    //! get input coupling regiontype
    virtual CouplingRegionType GetInputRegionType(UInt i)
    { return inputInterfaces_[i]->regionType; }

    //! get input coupling region nodes
    virtual void GetInputNodes(UInt i, StdVector<UInt>* &nodes)
    { nodes  = &(inputInterfaces_[i]->nodes);}

    //! get input coupling region elements
    virtual void GetInputElements(UInt i, StdVector<Elem *>*  &elements)
    { elements = &(inputInterfaces_[i]->elements);}

    //! get input neighbour region
    virtual void GetInputNeighbourRegion(UInt i, StdVector<std::string>* &regions)
    { regions = &(inputInterfaces_[i]->neighInputRegions);}

    //! get input coupling values
    virtual void GetInputValues(UInt i, SingleVector* &values)
    { values = (inputInterfaces_[i]->values);}

    //! get input coupling values
    virtual void GetInputOldValues(UInt i, SingleVector* &values)
    { values = (inputInterfaces_[i]->oldValues);}

    //! get input coupling values of step tn-1
    virtual void GetInputValues_tn_1(UInt i, SingleVector* &values_tn_1)
    { values_tn_1 = (inputInterfaces_[i]->values_tn_1);}

    //! get input coupling values of step tn-2
    virtual void GetInputValues_tn_2(UInt i, SingleVector* &values_tn_2)
    { values_tn_2 = (inputInterfaces_[i]->values_tn_2);}

    //! get input coupling values of step tn-3
    virtual void GetInputValues_tn_3(UInt i, SingleVector* &values_tn_3)
    { values_tn_3 = (inputInterfaces_[i]->values_tn_3);}

    //! get input coupling values dof
    virtual UInt GetInputDof(UInt i)
    { return inputInterfaces_[i]->dof; }

    //! get input coupling nodes size
    virtual UInt GetInputNumNodes(UInt i)
    { return inputInterfaces_[i]->numNodes; }

    //! get input coupling elems size
    virtual UInt GetInputNumElems(UInt i)
    { return inputInterfaces_[i]->numElems; }

    //! get input coupling norm type
    virtual NormType GetInputNormType(UInt i)
    { return inputInterfaces_[i]->normtype; }

    //! get input coupling epsilon
    virtual Double GetInputEpsilon(UInt i)
    { return inputInterfaces_[i]->epsilon; }

    // ------------- output coupling -----------
  
    //! get output coupling type
    virtual CouplingOutputType GetOutputType(UInt i)
    { return outputTypes_[i];}

    //! get output coupling quantity
    virtual SolutionType GetOutputQuantity(UInt i)
    { return outputQuantities_[i]; }

    //! get output coupling region
    virtual void GetOutputRegions(UInt i, StdVector<std::string> &regions)
    { regions = outputInterfaces_[i]->regions; }

    //! get output coupling regiontype
    virtual CouplingRegionType GetOutputRegionType(UInt i)
    { return outputInterfaces_[i]->regionType; }
 
    //! get output coupling region nodes
    virtual void GetOutputNodes(UInt i, StdVector<UInt>* &nodes)
    { nodes  = &(outputInterfaces_[i]->nodes);}

    //! get output coupling region elements
    virtual void GetOutputElements(UInt i, StdVector<Elem *>* &elements)
    { elements = &(outputInterfaces_[i]->elements);}
  
    //! get output neighbour elements
    //virtual void GetOutputNeighbourElems(UInt i, StdVector<Elem *>*  &elements)
    //{ elements = &(outputInterfaces_[i]->neighbours);}
  
    //! get output neighbour region
    virtual void GetOutputNeighbourRegion(UInt i, StdVector<std::string>* &regions)
    { regions = &(outputInterfaces_[i]->neighInputRegions);}
    
    //! get output coupling region materials
    virtual void GetMaterials(UInt i, StdVector<BaseMaterial *>* &mat)
    { mat = &(outputInterfaces_[i]->materials);} 
  
    //! get output coupling region materials of opposite pde
    virtual void GetOppositeMaterials(UInt i, StdVector<BaseMaterial *>* &mat)
    { mat = &(outputInterfaces_[i]->oppositePdeMaterials);}

    //! get output coupling values
    virtual void GetOutputValues(UInt i, SingleVector* &values)
    { values = (outputInterfaces_[i]->values);}

    //! get old output coupling values
    virtual void GetOutputOldValues(UInt i, SingleVector* &values)
    { values = (outputInterfaces_[i]->oldValues);}

    //! get output coupling values dof
    virtual UInt GetOutputDof(UInt i)
    { return outputInterfaces_[i]->dof; }

    //! get output coupling number nodes
    virtual UInt GetOutputNumNodes(UInt i)
    { return outputInterfaces_[i]->numNodes; }

    //! get output coupling number elems
    virtual UInt GetOutputNumElems(UInt i)
    { return outputInterfaces_[i]->numElems; }

    //! get output coupling norm type
    virtual NormType GetOutputNormType(UInt i)
    { return outputInterfaces_[i]->normtype; }

    //! get output coupling epsilon
    virtual Double GetOutputEpsilon(UInt i)
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


#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class PDECoupling
  //! 
  //! \purpose 
  //! 
  //! \collab 
  //! 
  //! \implement 
  //! 
  //! \status In use
  //! 
  //! \unused 
  //! 
  //! \improve
  //! 

#endif

} // end of namespace

#endif
