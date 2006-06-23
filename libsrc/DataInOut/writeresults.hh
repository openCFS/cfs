#ifndef FILE_WRITERESULTS_2002
#define FILE_WRITERESULTS_2002

#include "Domain/grid.hh"
#include "Utils/nodestoresol.hh"
#include "Utils/elemstoresol.hh"
#include "Utils/StdVector.hh"

// maximum number of open history files
#define MAX_NUM_HIST_FILES 1000

namespace CoupledField {

  //! Base class for writing results
  class WriteResults {

    //! Enum describing entities
    typedef enum {NO_ENTITIY,ELEMENT,NODE} entityType;
    
  public:
    //! Constructor
    WriteResults(const Char * const filename);

    //! initialization with grid
    //! \param ptgrid pointer to grid object
    virtual void Init(Grid * aptgrid)=0;

    //! Default destructor
    virtual ~WriteResults();

    //! write information about grid  in file
    virtual void WriteGrid()=0;
 
    // *************************
    //     TRANSIENT SECTION
    // *************************
  
    //! write element solution vector (transient/static)
    /*!
      \param data vector with data (ex. value of an error for the cell)
      \param step step of calculation
      \param time time of calculation
    */
    virtual void WriteNodeSolutionTransient(const NodeStoreSol<Double>& data, 
                                            const UInt step, 
                                            const Double time) = 0;

    //! write element solution vector (transient/static)
    /*!
      \param data vector with data (ex. value of an error for the cell)
      \param step step of calculation
      \param time time of calculation
    */
    virtual void WriteElemSolutionTransient(const ElemStoreSol<Double>& data, 
                                            const UInt step, 
                                            const Double time) = 0;

    //! write node history vector (transient/static)
    /*!
      \param data vector with data (ex. value of an error for the node)
      \param step step of calculation
      \param time time of calculation
    */
    virtual void WriteNodeHistoryTransient(const NodeStoreSol<Double>& data, 
                                           const UInt step, 
                                           const Double time);

    //! write element history vector (transient/static)
    /*!
      \param data vector with data (ex. value of an error for the cell)
      \param step step of calculation
      \param time time of calculation
    */
    virtual void WriteElemHistoryTransient(const ElemStoreSol<Double>& data, 
                                           const UInt step, 
                                           const Double time);
  

    // *************************
    //     HARMONIC SECTION
    // *************************

    //! write element solution vector (harmonic)
    /*!
      \param data vector with data (ex. value of an error for the cell)
      \param step step of calculation
      \param frequency frequency of exciting function
      \param format format for writing complex solution
      (real-imag/amplitude-phase)
    */
    virtual void WriteNodeSolutionHarmonic(const NodeStoreSol<Complex>& data, 
                                           const UInt step,
                                           const Double frequency,
                                           const ComplexFormat format) = 0;

    //! write element solution vector (harmonic)
    /*!
      \param data vector with data (ex. value of an error for the cell)
      \param step step of calculation
      \param frequency frequency of exciting function
      \param frequencyStep step of calculation
      \param format format for writing complex solution
      (real-imag/amplitude-phase)
    */
    virtual void WriteElemSolutionHarmonic(const ElemStoreSol<Complex>& data, 
                                           const UInt step,
                                           const Double frequency,
                                           const ComplexFormat format) = 0;

    //! write nodal history vector (harmonic)
    /*!
      \param data vector with data (ex. value of an error for the cell)
      \param step step of calculation
      \param frequency frequency of exciting function
      \param frequencyStep step of calculation
      \param format format for writing complex solution
      (real-imag/amplitude-phase)
    */
    void WriteNodeHistoryHarmonic(const NodeStoreSol<Complex>& data, 
                                  const UInt step,
                                  const Double frequency,
                                  const ComplexFormat format);

    //! write element history vector (harmonic)
    /*!
      \param data vector with data (ex. value of an error for the cell)
      \param step step of calculation
      \param frequency frequency of exciting function
      \param frequencyStep step of calculation
      \param format format for writing complex solution
      (real-imag/amplitude-phase)
    */
    virtual void WriteElemHistoryHarmonic(const ElemStoreSol<Complex>& data, 
                                          const UInt step,
                                          const Double frequency,
                                          const ComplexFormat format);
    

    //! to open new file for printing results only for GMV

    //! \param number number for output-file (ex. result.gmv001)
    virtual void OpenFile(const UInt number)
    { Error("Not implemented",__FILE__,__LINE__);}


    //! writes a matrix of the form ptCoordX ptCoorY ptCoordZ SolDof1 SolDof2 ...
    /*!
      \param nrDofs nr of degrees of freedom in solution
    */
    void WriteSolMatrix(Grid * ptgrid,
                        const Vector<Double> sol, 
                        const std::string matFileName, const UInt nrDofs=1);


  protected:

    //! Initializes the history data
    void InitHistoryData( StdVector<SolutionType> & quantities,
                          StdVector<StdVector<UInt> > & entitiesPerQuant,
                          StdVector<StdVector<std::ofstream*> > & histEntFiles,
                          entityType entType );
      
    //! Convertes enum SolutionType to string
    virtual std::string SolutionTypeToString(const SolutionType type) const {
      Error( "Not implemented here", __FILE__, __LINE__ );
      return "badluck";
    }

    //! name of file for output results
    std::string namefile_;

    //! pointer to Grid
    Grid * ptGrid_;

    // ************************************
    //  Section dealing with history files 
    // ************************************

    //! indicator: print history file or not
    bool NeedHistory_;

    //! Total number of open history files
    UInt totalNumHistFiles_;

    //@{ 
    // \name Data structures for nodal history
    //! vector with type of output values
    //! for which a history file is written
    StdVector<SolutionType> histNodeQuantities_;
  
    //! history nodes per output quantity
    StdVector<StdVector<UInt> > histNodesPerQuant_;

    //! pointer to ofstream with history information
    StdVector<StdVector<std::ofstream*> >  histNodeFiles_;
    //@}
    
    //@{
    // \name Data structures for element history
    //! vector with type of output values
    //! for which a history file is written
    StdVector<SolutionType> histElemQuantities_;
  
    //! history nodes per output quantity
    StdVector<StdVector<UInt> > histElemsPerQuant_;

    //! pointer to ofstream with history information
    StdVector<StdVector<std::ofstream*> >  histElemFiles_;

    //! intialize the management of history files
    virtual void InitHistoryFiles();

  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class WriteResults
  //! 
  //! \purpose
  //! This is the base class for writing simulation results. The methods of
  //! this class are purely virtual and must be over-written by the derived
  //! classes that handle the different types of output formats/files.
  //!
  //! \collab 
  //! The only object of this class gets instanciated in DefineFiles at the 
  //! bginning of a simulation run. The results themselves are passed to this
  //! class in the different SinglePDEs in the method WriteResultsInFile.
  //! 
  //! \implement 
  //! 
  //! \status In use
  //! 
  //! \unused 
  //! 
  //! \improve
  //! - Right now we open a new file for each history-node and quantity. 
  //! Since all of these files occupie a file descriptor during the whole
  //! simulation, there can occur strange errors when the operating systems
  //! runs out of file handlers (normally only 1024 are allowed per process!).
  //! Here  we need an appropriate check or we have to store all history 
  //! results in a vector and write them out only at the end of a simulation.

#endif

} // end of namespace

#endif
