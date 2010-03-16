// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_MPCCIPDE_NEW
#define FILE_MPCCIPDE_NEW

#include "SinglePDE.hh" 

#include <def_use_mpcci.hh>

#ifdef MpCCI
#include <MpCCIcpl/MpCCIexch.hh>
#endif

namespace CoupledField
{

  //! Class for coupling a pde via MpCCI

  class MpcciPDE : public SinglePDE {

  public:

    //! Constructor. here we read integration parameters
    /*!
      \param 
      \param aGrid pointer to grid
    */
    MpcciPDE( Grid* aptgrid, PtrParamNode paramNode );

    //! Deconstructor
    virtual ~MpcciPDE(){};

    virtual void Init(UInt sequenceStep = 0 );

    virtual void PreparePDE4Computation();

    //! define algebraic system 
    void DefineAlgSys();

    //! define all (bilinearform) integrators needed for this pde
    virtual void DefineIntegrators();

    //! define the SoltionStep-Driver
    virtual void DefineSolveStep();

    //! initialize time stepping: 
    //! nothing to do in mpcci!

    //! set current time step

   
    // ======================================================
    // POSTPROCESSING SECTION
    // ======================================================

    virtual void WriteGeneralPDEdefines() {};

  
    void WriteRestart(const UInt nstep, UInt totalUnknowns=0);


    void ReadRestart(UInt &startStep, UInt totalUnknowns=0);

    // ======================================================
    // COUPLING SECTION
    // ======================================================

    //! initalize PDE coupling
    virtual void InitCoupling(PDECoupling * Coupling);

    virtual void CalcInputCoupling();

    //! calculate coupling terms
    virtual void CalcOutputCoupling();

    //! returns if PDE can compute the quantity
    virtual bool HasOutput(SolutionType output);
  

  protected:

    //!Get an Integer array including all nodes of a subdomain.
    //!Thereby for each elements the connecting nodes are picked up and stored 
    //!sequentially in localNodes_. Thereby it is possible that one node
    //!is mentioned more then once. Because one node can belong to more elements.
    void GetNodesOfSubdomain();


    //!set up an bool Matrix with the mapping which pdenode belongs to which subdomain
    void SetupNodesSubdomainsMapping();


    StdVector<StdVector<Elem*> > F_Interface_; //!<vector of vectors conaining Elements with acting force
    StdVector<StdVector<StdVector<ShortInt> > > isBoundaryNode_; //!< vector containing flag array for element boundary nodes
    StdVector<StdVector<StdVector<UInt> > > elemNodeToCouplingNode_; //!< assigns each coupling element node the according Coupling Node number
    
    // *****************
    //  POSTPROCESSING
    // *****************

    // for check: own solver
    bool SolverCFS_; //<! parameter indicator: true, if you want to use Solver CFS. reading from config-file
    Matrix<Double> sysmat_;
    Vector<Double> vecrhs_;

  private:
    //!MpCCI
    StdVector<UInt> mapSD_;
    UInt * meshId_;
    UInt NumMeshIds_;
#ifdef MpCCI
    MpCCIexch * ptMpCCIexch_;
#endif
    bool flagFirstTimeStep_; // flag for first time stewp
    UInt MpCCInodes_; //<! number of FE-nodes for MpCCI-domain
    std::string MpCCIType_; // <! Coupling Type: shell:= pressure from bioth sides; solid:= pressure from one side

    // Matrix containig the local PDE node numbers of all subdomains of the PDE
    UInt ** localNodes_;
    UInt * numOfNodesInSD_;
    StdVector<StdVector<bool> > NodeBelongsToSD_;

  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class MpcciPDE
  //! 
  //! \purpose 
  //! This class is derived from class BasePDE. It is used for coupling CFS++ 
  //! via MpCCI to FASTEST-3D. The class is thereby a dummy class in which
  //! nothing needs to be solved here just data transfer is handled.
  //! Currently only a coupled computation with mechanics is possible. 
  //! Within the Method CalcInputCoupling the 
  //! mechanical displacements from CFS++ are delivered to FASTEST-3D.
  //! In the method CalcOutputCoupling the fluid dynamical forces due to
  //! pressure and shear are send to CFS++.
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
