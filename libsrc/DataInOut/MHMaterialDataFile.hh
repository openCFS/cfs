/**************************************************************************/
/* File:   MHMaterialDataFile.hh                                          */
/* Author: Tom Lahmer                                                     */
/* Date:   06.06.2005                                                     */ 
/*         Rev. 16.12.1999                                                */
/*                                                                        */
/* Creates, writes and reads exponent arrays for multiharmonic ansatz     */
/*                                                                        */
/**************************************************************************/

#include <General/environment.hh>
#include <Matrix/matrix.hh>


namespace CoupledField 
{

  
  class MHMaterialDataFile {

  public:
    //! Constructor
    MHMaterialDataFile();

    //! Destructor
    ~MHMaterialDataFile();

    std::ofstream * indexSetFile;
    std::ofstream * indexSetCountFile;

    std::ifstream * indexSetFileIn;
    std::ifstream * indexSetCountFileIn;

    void createFiles();
    //! Generates and writes an index set to file
    void computeIndexSet(Integer PP, UInt PPmax, UInt nrMultHarms_);

    //! Retrieves values from index file and provides array of exponents
    void getExponentArray(Matrix<Integer> & exponent, Matrix<Integer> & count, UInt pp, UInt p, Integer delta);
    
  protected:

  private:

  }; // end of class MHMaterialDataFile

#ifdef DOXYGEN_DETAILED_DOC
  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class MHMaterialDataFile
  //! 
  //! \purpose 
  //! Provides different types of material and parameters needed for the
  //! multiharmonic Driver. E.g. the exponents for the multinomial-theorem
  //! will be generated, stored into according files
  //! 
  //! \collab 
  //! 
  //! \implement 
  //! 
  //! \status Deprecated. 
  //! under construction
  //! 
  //! \unused 
  //! 
  //! \improve
  //! 

#endif

}



