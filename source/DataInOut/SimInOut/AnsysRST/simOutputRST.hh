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

  //! Class for writing results in the GiD postprocessing format
  class SimOutputRST: virtual public SimOutput {
    
  public:

    //! Constructor
    SimOutputRST( const std::string& fileName, ParamNode * outputNode );
  
    //! Destructor
    virtual ~SimOutputRST();
  
    //! Initialize class
    void Init( Grid* grid );

    //! Write grid definition in file
    void WriteGrid();

    //! Register result (within one multisequence step)
    void RegisterResult( shared_ptr<BaseResult> sol,
                         UInt saveBegin, UInt saveInc,
                         UInt saveEnd );
    
    //! Begin single analysis step
    void BeginStep( UInt stepNum, Double stepVal );

    //! Add result to current step
    void AddResult( shared_ptr<BaseResult> sol );

    //! End single analysis step
    void FinishStep( );

    void Finalize();

  private:

    //! Write nodal coordinates
    void WriteNodes();
    
    //! Write element declarations
    void WriteElements();
    
    //! Write transient results on nodes or elements
    void WriteNodeElemDataTrans( const Vector<Double> & var, 
                                 const StdVector<std::string> & dofNames,
                                 const std::string& name, 
                                 ResultInfo::EntryType entryType,
                                 ResultInfo::EntityUnknownType entityType,
                                 Double time );

    //! Write harmonic results on nodes or elements
    void WriteNodeElemDataHarm( const Vector<Complex>  & var, 
                                const StdVector<std::string> & dofNames,
                                std::string name, 
                                ResultInfo::EntryType entryType,
                                ResultInfo::EntityUnknownType entityType,
                                Double freq, 
                                ComplexFormat outputFormat );
    
    //! Reorder respectively degenerate element connectivity

    //! This method takes the connectivity of an element, and replaces it
    //! with the degenerated connectivity of the element if necessary.
    //! The method returns the ANSYS element type.
    Integer ReorderConnect4Ansys(FEType eType,
                                 StdVector<UInt>& connect,
                                 Integer* connectANSYS);

    enum ansys_elementtype
    {
      IELCSZ = 150,
      IETYP=1,
      JETYP=2,
      KYOP1=3,
      KYOP2=4,
      KYOP3=5,
      KYOP4=6,
      KYOP5=7,
      KYOP6=8,
      KYOP7=9,
      KYOP8=10,
      KYOP9=11,
      KYOP10=12,
      KYOP11=13,
      KYOP12=14,
      //  --- definitions of KYOP1 to KYOP12 must be continuous and in
      // ---- ascending order as other parts of code presume this

      //*** maximum number of elements in ansys library (largest JETYP)
      MAXLIB = 200,

      // define output indices

      KDIM=21,
      ISHAP=22,
      IDEGEN=23,
      JSHELL=24,
      MNODE=25,
      MSHLIM=26,
      KCOIN=27,
      KELSTO=28,
      KLAYER=29,
      JPIPE=30,
      NMLYMX=31,
      SHCHEK=32,
      //      ---- reshuffle and rename for layers.
      KDOFS=34,
      NMVCTF=35,
      NOTRAN=36,
      KDOFF=37,
      NMDRLC=39,
      NMTRLC=40,
      KNORLC=41,
      LCANGL=42,
      KXYZRC=43,
      MATRQD=47,
      MATISO=48,
      MATMUL=49,
      ISTRES=50,
      TBRQD=51,
      NCOMPN=52,
      KSUPER=55,
      MTSYM=56,
      POSDF=57,
      INSS=58,
      INLD=59,
      KSEKU=60,
      NMNDMX=61,        // maximum number of nodes per element
      NMNDMN=62,
      NMNDST=63,
      NMDFPN=64,
      NMNDAC=65,
      MATRXS=66,
      NMNDNE=67,
      KTANS =68,
      JMASS=70,
      JSECT=71,
      JORIENT=72,
      NMPTSF=74,
      NMPRES=75,
      NMCONV=76,
      NMIMSF=77,
      NMMFLU=78,
      NMFGSF=79,
      NMCDSF=80,
      NMRAD=81, 
      NMPORT=82,
      NMTEMP=83,
      NMFLNC=84,
      NMHTGN=85,
      NMMSCN=86,
      NMVTDP=87,
      NMCRDN=88,
      NMVLTG=89,
      NMVLCD=90,
      NMFRSF=91,
      NMRDSF=92,
      NMVFRC=93,
      NMNDNO=94,  // number of corner nodes
      KCONIT=95,
      KLINTP=96,
      NMLINX=97,
      NMLINR=98,
      NMSHLR=99,
      MEASTR=100,
      KINITS=101,
      NMBFEF=102,
      NMBFH=103,
      NMBFPRT=104,
      KRDB=105,
      KNORM=106,
      NMSFST=108,
      NMSMIS=109,
      NMNMIS=110,
      NMNMUP=111,
      NCPTM1=112,
      NCPTM2=113,
      NFILER=114,
      KMAGC=115,
      EDGEEL=116,
      // -- at next restructure, put NMSMIS,NMNMIS,and NMNMUP at
      // -- 108,109,110 to facilitate study of echtst output - pck
      JSTPR=117,
      NMICVF=118,
      NMINND=119, 
      NMSHLD=120, 
      NMPERB=121,
      NMSSVR=122,
      NMFSVR=123,
      NMMSVR=124,
      NBSVSP=125,
      NMNSVR=126,
      NMPSVR=127,
      NMCSVR=128,
      NMUSVR=129,
      NMOSVR=131,
      ELFORCE=132,   
      KSWTEI=133,
      //  ---- at next restructure, put KSWTEI and KSWTTS together - pck
      NAMRF1=134,
      NAMRF2=135, 
      ADADES=137,
      NOESAV=138,   
      NLSSAD=139,
      JDEAD=140,
      KSWTTS=141,
      KPSO=142,
      JNONLR=143,
      JSTAT=144,
      NOEMAT=145,
      KFLUID=146,
      KSURF=147,
      KTRG=148,
      PELEM=149,    
      NOELEM=150
    };

    enum ansys_nodal_dof
    {
      UX   =  1,
      UY   =  2,
      UZ   =  3,
      ROTX =  4,
      ROTY =  5,
      ROTZ =  6,
      AX   =  7,
      AY   =  8,
      AZ   =  9,
      VX   = 10,
      VY   = 11,
      VZ   = 12,

      PRES = 19,
      TEMP = 20,
      VOLT = 21,
      MAG  = 22,
      ENKE = 23,
      ENDS = 24,
      EMF  = 25,
      CURR = 26,
      SP01 = 27,
      SP02 = 28,
      SP03 = 29,
      SP04 = 30,
      SP05 = 31,
      SP06 = 32
    };

    typedef struct{
      Integer elementtypid;
      Integer pos;
      Integer ielc[IELCSZ];
    } AnsysElementType;

    void SetElementType42(AnsysElementType& eltype,
                          Integer TypeNumber);
    void SetElementType45(AnsysElementType& eltype,
                          Integer TypeNumber);
    void SetElementType82(AnsysElementType& eltype,
                          Integer TypeNumber);
    void SetElementType95(AnsysElementType& eltype,
                          Integer TypeNumber);

    
    //! Map with result objects for each result type
    ResultMapType resultMap_;

    //! Check if class is initialized
    bool isInitialized_;

    Double ansysNode_[6];
    Integer connectANSYS_[27];
    Integer elemData_[10];
    std::map<Integer, Integer> ansysType2TypeId_;
    std::map<Integer, Integer> ansysType2NumNodes_;


  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class SimOutputRST
  //! 
  //! \purpose 
  //! This class provides the interface for writing meshes and elements to the
  //! GiD own postprocessing files .post.msh and .post.res.
  //! This format is capable of visualizing volume elements, surface elements
  //! and named nodes. 
  //! /Note GiD has some restrictions regarding the visualization of element
  //! resuls on different element types. Therefore, if different types of 
  //! elements are contained in the grid (e.g. Quad/Tria, Tet/Hex/Pyra), all
  //! elements are treated as degenerated quadrilaterals or hexahedras. 
  //! This might cause some distortions in the resulting visualization, 
  //! especially if quadratic elements are used.
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
