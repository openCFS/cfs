// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_SIM_OUTRST_HH
#define FILE_CFS_SIM_OUTRST_HH

#include <Domain/grid.hh>
#include <DataInOut/simOutput.hh>
#include <Domain/resultInfo.hh>

namespace CoupledField
{

  //! Class for writing results in ANSYS RST postprocessing format
  class SimOutputRST: virtual public SimOutput {
    
  public:

    // We use 2D/3D ANSYS mechanic elements.
    enum UsedAnsysElementTypes
    {
      PLANE42 = 42,
      SOLID45 = 45,
      PLANE82 = 82,
      SOLID95 = 95
    };

    // These are the 150 members of the ANSYS element header
    enum AnsysElementTypeDescr
    {
      IELCSZ = 150,
      IETYP=1,    JETYP=2,    KYOP1=3,    KYOP2=4,    KYOP3=5,
      KYOP4=6,    KYOP5=7,    KYOP6=8,    KYOP7=9,    KYOP8=10,
      KYOP9=11,   KYOP10=12,  KYOP11=13,  KYOP12=14,  
      //  --- definitions of KYOP1 to KYOP12 must be continuous and in
      // ---- ascending order as other parts of code presume this

      //*** maximum number of elements in ansys library (largest JETYP)
      MAXLIB = 200,

      // define output indices
      KDIM=21,    ISHAP=22,   IDEGEN=23,  JSHELL=24,  MNODE=25,
      MSHLIM=26,  KCOIN=27,   KELSTO=28,  KLAYER=29,  JPIPE=30,
      NMLYMX=31,  SHCHEK=32,  
      //      ---- reshuffle and rename for layers.
      KDOFS=34,   NMVCTF=35,  NOTRAN=36,  KDOFF=37,   NMDRLC=39,
      NMTRLC=40,  KNORLC=41,  LCANGL=42,  KXYZRC=43,  MATRQD=47,
      MATISO=48,  MATMUL=49,  ISTRES=50,  TBRQD=51,   NCOMPN=52,
      KSUPER=55,  MTSYM=56,   POSDF=57,   INSS=58,    INLD=59,
      KSEKU=60,   NMNDMX=61, // maximum number of nodes per element
      NMNDMN=62,  NMNDST=63,  NMDFPN=64,  NMNDAC=65,  MATRXS=66,
      NMNDNE=67,  KTANS =68,  JMASS=70,   JSECT=71,   JORIENT=72,
      NMPTSF=74,  NMPRES=75,  NMCONV=76,  NMIMSF=77,  NMMFLU=78,
      NMFGSF=79,  NMCDSF=80,  NMRAD=81,   NMPORT=82,  NMTEMP=83,
      NMFLNC=84,  NMHTGN=85,  NMMSCN=86,  NMVTDP=87,  NMCRDN=88,
      NMVLTG=89,  NMVLCD=90,  NMFRSF=91,  NMRDSF=92,  NMVFRC=93,
      NMNDNO=94,  // number of corner nodes
      KCONIT=95,  KLINTP=96,  NMLINX=97,  NMLINR=98,  NMSHLR=99,
      MEASTR=100, KINITS=101, NMBFEF=102, NMBFH=103,  NMBFPRT=104,
      KRDB=105,   KNORM=106,  NMSFST=108, NMSMIS=109, NMNMIS=110,
      NMNMUP=111, NCPTM1=112, NCPTM2=113, NFILER=114, KMAGC=115,
      EDGEEL=116, 
      // -- at next restructure, put NMSMIS,NMNMIS,and NMNMUP at
      // -- 108,109,110 to facilitate study of echtst output - pck
      JSTPR=117,  NMICVF=118, NMINND=119, NMSHLD=120, NMPERB=121,
      NMSSVR=122, NMFSVR=123, NMMSVR=124, NBSVSP=125, NMNSVR=126,
      NMPSVR=127, NMCSVR=128, NMUSVR=129, NMOSVR=131, ELFORCE=132,   
      KSWTEI=133, 
      //  ---- at next restructure, put KSWTEI and KSWTTS together - pck
      NAMRF1=134, NAMRF2=135, ADADES=137, NOESAV=138, NLSSAD=139,
      JDEAD=140,  KSWTTS=141, KPSO=142,   JNONLR=143, JSTAT=144,
      NOEMAT=145, KFLUID=146, KSURF=147,  KTRG=148,   PELEM=149,    
      NOELEM=150
    };

    // ANSYS nodal DOFs
    enum AnsysNodalDof
    {
      UX   =  1, UY   =  2, UZ   =  3,
      ROTX =  4, ROTY =  5, ROTZ =  6,
      AX   =  7, AY   =  8, AZ   =  9,
      VX   = 10, VY   = 11, VZ   = 12,

      PRES = 19,
      TEMP = 20,
      VOLT = 21,
      MAG  = 22,
      ENKE = 23,
      ENDS = 24,
      EMF  = 25,
      CURR = 26,
      SP01 = 27, SP02 = 28, SP03 = 29,
      SP04 = 30, SP05 = 31, SP06 = 32,
      UNKN = -1
    };

    // ANSYS element DOFs
    enum AnsysElemDof
    {
      EMS= 1, ENF= 2, ENS= 3, ENG= 4, EGR= 5,
      EEL= 6, EPL= 7, ECR= 8, ETH= 9, EUL=10,
      EFX=11, ELF=12, EMN=13, ECD=14, ENL=15,
      EHC=16, EPT=17, ESF=18, UNK1=19, ETB=20,
      ECT=21, EXY=22, EBA=23, ESV=24, UNK2=25,
      UNKE=-1
    };
    
    // Structure which can hold local element type and header. 
    typedef struct{
      Integer elementtypid;
      Integer ielc[IELCSZ];
    } AnsysElementType;

  public:

    //! Constructor
    SimOutputRST( const std::string& fileName, ParamNode * outputNode );
  
    //! Destructor
    virtual ~SimOutputRST();
  
    //! Initialize class
    void Init( Grid* grid, bool printGridOnly );

    //! Write grid definition in file
    void WriteGrid();

    //! Begin multisequence step
    virtual void BeginMultiSequenceStep( UInt step,
                                         AnalysisType type,
                                         UInt numSteps);

    //! Register result (within one multisequence step)
    void RegisterResult( shared_ptr<BaseResult> sol,
                         UInt saveBegin, UInt saveInc,
                         UInt saveEnd,
                         bool isHistory );
    
    //! Begin single analysis step
    void BeginStep( UInt stepNum, Double stepVal );

    //! Add result to current step
    void AddResult( shared_ptr<BaseResult> sol );

    //! End single analysis step
    void FinishStep( );

    //! End multisequence step
    virtual void FinishMultiSequenceStep( );

    //! Deallocate last resources
    void Finalize();

  private:

    //! Write nodal coordinates
    void WriteNodes();
    
    //! Write element declarations
    void WriteElements();
    
    //! Prepare transient nodal data for passing to ANSYS Fortran 90 API
    void PrepareNodeData( const Vector<Double> & var,
                        std::vector<Double> & nodeValuesReal,
                        SolutionType solType,
                        UInt numDOFs,
                        const ResultInfo::EntryType entryType);

    //! Prepare transient element data for passing to ANSYS Fortran 90 API
    void PrepareElemData( const Vector<Double>  & var, 
                          SolutionType solType,
                          UInt numDOFs,
                          ResultInfo::EntryType entryType);


    //! Prepare complex nodal data for passing to ANSYS Fortran 90 API
    void PrepareNodeData( const Vector<Complex> & var,
                          std::vector<Double> & nodeValuesReal,
                          std::vector<Double> & nodeValuesImag,
                          SolutionType solType,
                          UInt numDOFs,
                          const ResultInfo::EntryType entryType);
    
    //! Prepare complex element data for passing to ANSYS Fortran 90 API
    void PrepareElemData( const Vector<Complex>  & var, 
                          SolutionType solType,
                          UInt numDOFs,
                          ResultInfo::EntryType entryType);
    
    //! Reorder respectively degenerate element connectivity

    //! This method takes the connectivity of an element, and replaces it
    //! with the degenerated connectivity of the element if necessary.
    //! The method returns the ANSYS element type.
    Integer ReorderConnect4Ansys(FEType eType,
                                 StdVector<UInt>& connect,
                                 Integer* connectANSYS);

    //! Map CFS/NACS nodal DOFs to ANSYS DOFs
    void MapInternal2ANSYSNodeDof(SolutionType solType);

    //! This map stores the mapping from CFS/NACS nodal DOFs to ANSYS DOFs.
    std::map <SolutionType, AnsysNodalDof> internal2AnsysNodeDofMap_;

    //! This map stores the position in the output array of a given ANSYS DOF.
    std::map <AnsysNodalDof, UInt> ansysNodeDof2Idx_;

    //! Map CFS/NACS element DOFs to ANSYS DOFs
    void MapInternal2ANSYSElemDof(SolutionType solType);

    //! This map stores the mapping from CFS/NACS element DOFs to ANSYS DOFs.
    std::map <SolutionType, AnsysElemDof> internal2AnsysElemDofMap_;

    //! Initialize element header for PLANE42 element.
    void SetElementType42(AnsysElementType& eltype,
                          Integer TypeNumber);

    //! Initialize element header for SOLID45 element.
    void SetElementType45(AnsysElementType& eltype,
                          Integer TypeNumber);

    //! Initialize element header for PLANE82 element.
    void SetElementType82(AnsysElementType& eltype,
                          Integer TypeNumber);

    //! Initialize element header for SOLID95 element.
    void SetElementType95(AnsysElementType& eltype,
                          Integer TypeNumber);

    //! Map with result objects for each result type
    ResultMapType resultMap_;

    //! Check if the grid is to be printed.
    bool printGrid_;

    //! Due to the fact, that we get to know the number of Dofs
    //! only after all results are registered we set this flag to
    //! tell the WriteGrid method to print the grid now.
    bool printNow_;

    // Temporary array for ANSYS node coordinate
    Double ansysNode_[6];

    // Temporary array for ANSYS element connectivity
    Integer connectANSYS_[27];

    // Temporary array for ANSYS element data (material number, etc.)
    Integer elemData_[10];

    //! Map ANSYS element type to type id in RST file.
    std::map<Integer, Integer> ansysType2TypeId_;

    //! Map ANSYS element type number of corners.
    std::map<Integer, Integer> ansysType2NumNodes_;

    //! Map internal element types to ANSYS element types.
    std::map<FEType, Integer> elemType2AnsysType_;

    //! Map ANSYS element DOFs to their number of DOFs
    //! e.g. stress has 11 DOFs, etc.
    std::map<Integer, Integer> numAnsysElemDOFs_;

    //! Vector of ANSYS nodal DOF labels.
    std::vector<std::string> dofLabels_;

    //! Count number of nodal DOFs of simulation.
    Integer numNodeDOFs_;

    // Number of current time/frequency step.
    Integer stepNum_;

    // Time/Frequency of current step.
    Double stepVal_;
    
    // Offset for step number in case of multisequence analysis
    Integer stepNumOffset_;
    
    // Offset for step value in case of multisequence analysis
    Double stepValOffset_;

    // Number of time/frequency steps per multistep.
    Integer numSteps_;

    // Analysis type.
    AnalysisType analysisType_;

    // Multistep number.
    UInt msStep_;

    // Does analysis type provide complex results?
    bool complexResults_;

    // "time" for transient/static, "frequency" for harmonic/eigenfrequency.
    std::string timeFreq_;

    // Map a ANSYS element DOF to a temporary vector with real element values.
    std::map<AnsysElemDof, std::vector<Double> > tmpElemResultVecsReal_;

    // Map a ANSYS element DOF to a temporary vector with imag element values.
    std::map<AnsysElemDof, std::vector<Double> > tmpElemResultVecsImag_;
  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class SimOutputRST
  //! 
  //! \purpose 
  //! This class provides the interface for writing meshes and elements to
  //! ANSYS postprocessing files .rst. All element types are mapped to four
  //! different ANSYS structural analysis element types PLANE42, SOLID45,
  //! PLANE82 and SOLID95. These are linear/quadratic quads/hexas.
  //! All other element types become degenerated versions of the above types.
  //! At the moment this format is capable of visualizing volume elements, and
  //! surface elements. 
  //! The node solutions get mapped to ANSYS Dofs. There are 32 such Dofs
  //! available (see AnsysNodalDof) and each of them gets written for every
  //! node in the grid in each time step.
  //! This is a tremendous waste of memory and should be redone.
  //! Element results get mapped to one of the 25 ANSYS element result types
  //! (see. AnsysElemDof)
  //! /Note In the moment only transient results will be written. Trying
  //! to write harmonic results will cause an EXCEPTION!
  //! For details about ANSYS binary files see Chapter 1: Format of Binary data
  //! files in the Guide to Interfacing with ANSYS (Programmer's Manual for ANSYS)
  //! The routines used to actually write the files are described in Chapter 2:
  //! Accessing Binary data files. The contents of the files can either be dumped
  //! with the bintst utility compiled with CFS/NACS or in ANSYS with
  //! File->List->Binary Files...
  //! The source code of the unv2rst utility was a good source of info, too.
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
