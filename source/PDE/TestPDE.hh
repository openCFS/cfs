// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

// ===========================================================================
/*!
 *       \file     TestPDE.hh
 *       \brief    This PDE will be used to train new developers
 *
 *       \date     July 2, 2013
 *       \author   Manfred Kaltenbacher, TU Wien
 */
//============================================================================

#ifndef TESTPDE_HH_
#define TESTPDE_HH_

#include "SinglePDE.hh"
 

namespace CoupledField {

  // forward class declaration
  class BaseResult;
  class ResultHandler;
  class linElecInt;
  class LinearFormContext;
  class BaseBDBInt;
  class ResultFunctor;

  //! Class for TEST PDE
  class TestPDE: public SinglePDE
  {

  public:

    //!  Constructor. Here we read integration parameters.
    /*!
      \param aGrid pointer to grid
      \param paramNode pointer to the corresponding parameter node
      \param infoNode pointer for writing to info (*.info.xml) file
      \param simState pointer to state for stopping and restarting simulations
      \param domain pointer do domain object
    */
    TestPDE( Grid* aGrid, PtrParamNode paramNode,
             PtrParamNode infoNode,
             shared_ptr<SimState> simState, Domain* domain  );

    //! Destructor
    virtual ~TestPDE(){};

  protected:
    
    //! Initialize NonLinearities
    virtual void InitNonLin();

     //! \copydoc SinglePDE::CreateFeSpaces
    virtual std::map<SolutionType, shared_ptr<FeSpace> > 
    CreateFeSpaces( const std::string&  formulation,
                    PtrParamNode infoNode );

    //! define all (bilinearform) integrators needed for this pde
    void DefineIntegrators();

    //! define surface integrators needed for this pde
    void DefineSurfaceIntegrators( ){};

    //! Defines the integrators needed for ncInterfaces
    void DefineNcIntegrators() {};

    //! Define all RHS linearforms for load / excitation 
    void DefineRhsLoadIntegrators();
    
    //! define the SoltionStep-Driver
    void DefineSolveStep();

    //!  Define available postprocessing results
    void DefinePrimaryResults();

    //! Define available postprocessing results
    void DefinePostProcResults();

    //! Init the time stepping
    void InitTimeStepping();
    

  private:

    //! Obtain information on desired output quantities from parameter file

    //! This method is used to query the parameter handling object for the
    //! desired output quantities and translate their literal description into
    //! the internal format by setting the corresponding class attributes.
    //! The output quantities currently supported by the heat conduction PDE 
    //! are given in the following table. Here 'Keyword' and 'Result Type'
    //! refer to the XML parameter file, while 'Class Attribute' refers to the
    //! internal attribute of the TestPDE class that is set, if the
    //! keyword is specified.\n
    //! \todo Specification of ReadStoreResults for TestPDE!!!
    //void ReadStoreResults();

    
  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class TestPDE
  //! 
  //! \purpose 
  //! This class is derived from class SinglePDE.
  //! It is used for solving the test PDE on one time step.  
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
