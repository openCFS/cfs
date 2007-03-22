// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef ILDLK_FACTORISER_HH
#define ILDLK_FACTORISER_HH

#include "utils/utils.hh"
#include "matvec/matvec.hh"
#include "precond/ILDLPrecond/baseildlfactoriser.hh"


namespace OLAS {

  //! This class implements an incomplete LDL(k) factorisation

  //! This class implements an incomplete LDL factorisation. In this special
  //! variant the dropping of entries is based on their fill-level, which can
  //! be understood as their distance from the original matrix pattern.
  //!
  //! An %ILDLKFactoriser object reads the following parameters from its
  //! myParams_ object:
  //! \n\n
  //! <center>
  //!   <table border="1" width="80%" cellpadding="10">
  //!     <tr>
  //!       <td colspan="4" align="center">
  //!         <b>Parameters for %ILDLKFactoriser object</b>
  //!       </td>
  //!     </tr>
  //!     <tr>
  //!       <td align="center"><b>ID string in OLAS_Params</b></td>
  //!       <td align="center"><b>range</b></td>
  //!       <td align="center"><b>default value</b></td>
  //!       <td align="center"><b>description</b></td>
  //!     </tr>
  //!     <tr>
  //!       <td>ILDLPRECOND_level</td>
  //!       <td align="center">positive integer</td>
  //!       <td align="center">1</td>
  //!       <td>This integer value gives the maximal allowed fill-level for
  //!           entries in the factor \f$U=L^T\f$. All entries with a
  //!           fill-level that exceeds 'level' are discarded.
  //!       </td>
  //!     </tr>
  //!     <tr>
  //!       <td>ILDLPRECOND_saveLevels</td>
  //!       <td align="center">true or false</td>
  //!       <td align="center">false</td>
  //!       <td>If the value of this boolean parameter is true, then the
  //!           fill-level information computed for the matrix factor
  //!           \f$U=L^T\f$ is exported to a file in MatrixMarket format.
  //!       </td>
  //!     </tr>
  //!     <tr>
  //!       <td>ILDLPRECOND_levelsFileName</td>
  //!       <td align="center">string</td>
  //!       <td align="center">-</td>
  //!       <td>File name for exporting the fill-level information of the
  //!           factorisation matrix. This value is only required, when
  //!           ILDLPRECOND_saveLevels is 'true'.
  //!       </td>
  //!     </tr>
  //!   </table>
  //! </center>
  template <class T>
  class ILDLKFactoriser : public BaseILDLFactoriser<T> {

  public:

    //! Default constructor
    ILDLKFactoriser();

    //! Standard constructor
    ILDLKFactoriser( OLAS_Params *myParams, OLAS_Report *myReport = NULL );

    //! Default destructor
    ~ILDLKFactoriser();

    //! Incomplete \f$LDL^T\f$ factorisation of a square matrix

    //! This is the central method of the class. It computes an
    //! incomplete \f$LDL^T\f$ factorisation of a square matrix using the
    //! ILDL(k) variant as special dropping strategy.
    //! \param sysMat  the matrix to be factorised
    //! \param dataD   vector with entries of diagonal factor
    //! \param rptrU   row pointer array for U factor
    //! \param cidxU   column index array for U factor
    //! \param dataU   array containing entries of U factor
    //! \param newPattern a 'true' indicates, that the matrix pattern has
    //!                   changed, i.e. the factorisation will require a new
    //!                   pattern analysis phase
    void Factorise( SCRS_Matrix<T> &sysMat, std::vector<T> &dataD,
                    std::vector<UInt> &rptrU, std::vector<UInt> &cidxU,
                    std::vector<T> &dataU, bool newPattern );

  private:

    //! Level ...
    UInt level_;

    //! Perform analytic factorisation
    void FactoriseAnalytic( SCRS_Matrix<T> &sysMat, std::vector<T> &dataD,
                            std::vector<UInt> &rptrU, std::vector<UInt> &cidxU,
                            std::vector<T> &dataU );

    //! Perform numerical factorisation
    void FactoriseNumerics( SCRS_Matrix<T> &sysMat, std::vector<T> &dataD,
                            std::vector<UInt> &rptrU, std::vector<UInt> &cidxU,
                            std::vector<T> &dataU );

    //! Estimate the memory required for storing the factor matrix
    UInt EstimateFactorMemory( SCRS_Matrix<T> &sysMat, UInt *auxVec );

    //! Export fill-levels of entries in factor \f$U=L^T\f$
    void ExportFillLevels( const char *fname, std::vector<UInt> rptrU,
                           std::vector<UInt> cidxU, std::vector<UInt> fillU );


    //! Drop entries from factor row depending on their fill level

    //! This method turns the full LDL factorisation into an incomplete one by
    //! dropping entries from the factor rows. An entry is discarded when its
    //! fill level is larger than the allowed threshold given by the
    //! ILDLPRECOND_level steering parameter
    void DropEntries( UInt *listIDX, T *listLVAL, UInt *listLVL );

  };

}

#endif
