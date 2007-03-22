// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef ILDL0_FACTORISER_HH
#define ILDL0_FACTORISER_HH

#include "utils/utils.hh"
#include "matvec/matvec.hh"
#include "precond/ILDLPrecond/baseildlfactoriser.hh"


namespace OLAS {

  //! This class implements an incomplete LDL factorisation with no fill in

  template <class T>
  class ILDL0Factoriser : public BaseILDLFactoriser<T> {

  public:

    //! Default constructor
    ILDL0Factoriser();

    //! Standard constructor
    ILDL0Factoriser( OLAS_Params *myParams, OLAS_Report *myReport = NULL );

    //! Default destructor
    ~ILDL0Factoriser();

    //! Incomplete \f$LDL^T\f$ factorisation of a square matrix

    //! This is the central method of the class. It computes an
    //! incomplete \f$LDL^T\f$ factorisation of a square matrix according to
    //! A = (L+I)*D*(L+I)'.
    //! \param sysMat  the matrix to be factorised
    //! \param dataD   vector with entries of diagonal factor
    //! \param rptrU   row pointer array for U factor
    //! \param cidxU   column index array for U factor
    //! \param dataU   array containing entries of U factor
    //! \param newPattern always false, since matrix pattern remains
    void Factorise( SCRS_Matrix<T> &sysMat, std::vector<T> &dataD,
                    std::vector<UInt> &rptrU, std::vector<UInt> &cidxU,
                    std::vector<T> &dataU, bool newPattern );

  private:

    //! Perform numerical factorisation
    void FactoriseNumerics( SCRS_Matrix<T> &sysMat, std::vector<T> &dataD,
                            std::vector<UInt> &rptrU, std::vector<UInt> &cidxU,
                            std::vector<T> &dataU );

  };

}

#endif
