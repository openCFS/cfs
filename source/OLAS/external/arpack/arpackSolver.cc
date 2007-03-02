#include "arpackSolver.hh"

#include <limits>
#include "solver/solver.hh"
#include "precond/precond.hh"
#include "matvec/matvec.hh"

namespace OLAS {

  ArpackSolver::ArpackSolver() {
      ENTER_FCN( "ArpackSolver::ArpackSolver" );

      shiftAndInvert_ = TRUE;
      freqShift_ = 0.0;
      tolerance_ = 1.e-8;
      maxIterations_ = 1000;
      numArnoldiVec_ = 0;
      interface_ = NULL;
      which_ = "LM";

  }
  
  ArpackSolver::~ArpackSolver() {
      ENTER_FCN( "ArpackSolver::~ArpackSolver" );
    
      delete[] eigenValues_;
      delete[] eigenTolerances_;
      delete[] eigenVectors_;
      interface_ = NULL;

  }

  void ArpackSolver::Setup( ArpackMatInterface *matInterface, UInt size, 
              UInt numFreq, Double freqShift, Char* which, bool shiftMode ) {
      ENTER_FCN( "ArpackSolver::Setup" );

      size_ = size;

      // new reference to matrix interface
      interface_ = matInterface;
    
      // Save frequency parameters
      numFreq_ = numFreq;
      freqShift_ = freqShift;
      shiftAndInvert_ = shiftMode;
      which_ = which;

      // set default value for Arnoldi vectors
      numArnoldiVec_ = numFreq_*2;

      // eigenvalues, tolerances and vectors
      eigenValues_ = New Double [numArnoldiVec_];
      eigenTolerances_ = New Double [numArnoldiVec_];
      eigenVectors_ = New Double [numFreq_*size_];

      DebugOff();

  }

  void ArpackSolver::SetTolerance(Double tol) {
      //! set ARPACK convergence tolerance in case of user specified value
      tolerance_ = tol;
  }

  void ArpackSolver::SetIterations(Integer maxIt) {
      //! set ARPACK convergence tolerance in case of user specified value
      maxIterations_ = maxIt;
  }

  void ArpackSolver::SetNumVectors(Integer numVec) {
      //! set ARPACK convergence tolerance in case of user specified value
      numArnoldiVec_ = numVec;
  }

  UInt ArpackSolver::GetNev( ) {
      return numFreq_;
  }

  Char* ArpackSolver::GetWhich( ) {
      return which_;
  }

  Double ArpackSolver::GetShift( ) {
      return freqShift_;
  }

  Double ArpackSolver::GetTol( ) {
      return tolerance_;
  }

  UInt ArpackSolver::GetMaxit( ) {
      return maxIterations_;
  }

  UInt ArpackSolver::GetNcv( ) {
      return numArnoldiVec_;
  }

  UInt ArpackSolver::FindEigenvalues( ) {
      ENTER_FCN( "ArpackSolver::FindEigenvalues" );

      bool converged = FALSE;
      Integer  ido = 0, info = 0;

      // temp vector to store B*x
      Double *tempV = New Double [size_];
      // pointers to the position of x and y in workD
      Double *vecX, *vecY;

      // temp working space required in dsaupd
      Double *workD = New Double [3*size_];
      Double *residual = New Double [size_];
      Double *workL = New Double [numArnoldiVec_*(numArnoldiVec_+8)];
      Double *matrixV = New Double [size_*numArnoldiVec_];

      InitTempSpace(tempV, residual, workD, workL, matrixV);

      Integer lenWorkL = numArnoldiVec_*(numArnoldiVec_+8);

      Integer iparams[21], ipntr[21];
      for (Integer i=0; i<21; i++) {
           iparams[i] = 0;
           ipntr[i] = 0;
      }

      iparams[0] = 1;
      iparams[2] = maxIterations_;
      iparams[6] = 3;

      UInt itNum;
      for (itNum=1; itNum<=maxIterations_; itNum++) {

          ARPACK_DSAUPD(&ido, "G", (Integer*) &size_, which_, (Integer*) &numFreq_, 
                &tolerance_, residual, (Integer*) &numArnoldiVec_, matrixV, 
                (Integer*) &size_, iparams, ipntr,
                workD, workL, &lenWorkL, &info);

          switch (ido)  {

          case -1:
              // NOTE: addresses in workd are given as FORTRAN addreses
              // solve (A-shift*B)*y = B*x
              // x is returned by dsaupd in workd(ipntr[1]:ipntr[1]+size-1)
              // y is expected by dsaupd in workd(ipntr[2]:ipntr[2]+size-1)

             vecX = workD + --ipntr[0];
             vecY = workD + --ipntr[1];
             interface_->MultBV(vecX, tempV); 
             interface_->MultShiftOpV(tempV, vecY); 
             break;

          case 1:
              // solve (A-shift*B)*y = B*x
              // B*x is returned by dsaupd in workd(ipntr[3]:ipntr[3]+size-1)
              // y is expected by dsaupd in workd(ipntr[2]:ipntr(2]+size-1)
             vecX = workD + --ipntr[2];
             vecY = workD + --ipntr[1];
             interface_->MultShiftOpV(vecX,vecY);
             break;

          case 2:
              // calculate y = B*x
              // x is returned by dsaupd in workd(ipntr[1]:ipntr[1]+size-1)
              // y is expected by dsaupd in workd(ipntr[2]:ipntr[2]+size-1)
              vecX = workD + --ipntr[0];
              vecY = workD + --ipntr[1];
              interface_->MultBV(vecX,vecY);
              break;

          case 3:
              // algorithm requires shifts - > not implemented
              (*error) << "User required shifts not yet implemented!\n"
                       << "ARPACK requested ido=3 in DSAUPD!";
              Error( __FILE__, __LINE__ );
              break;

          case 99:
              converged = TRUE;
              break;

          default:
              // something went wrong -> arpack requires undefind ido
              (*error) << "ARPACK requested not implemented ido=" << ido 
                       << " in DSAUPD!";
              Error( __FILE__, __LINE__ );
              break;

          }

          // if convergence has been achieved, exit loop
          if (converged)
              break;

      }

      // check whether everything went well
      if (info<0) {
          (*error) << "Error reported in ritz value calculation!\n"
                   << "Errorcode reported is " << info << "!";
          Error( __FILE__, __LINE__ );
      }

      if ( itNum>maxIterations_ && info != 99 ) {
          (*error) << "No convergence achieved within maximum number of "
                   << "allowed iterations!\n"
                   << "Increase maxIt (=" << maxIterations_ << ") in xml file!";
          Error( __FILE__, __LINE__ );
      }

      // proceed with calculation of eigenvectors, values, and tolerances

      // NOTE: in the FORTRAN routine dseupd d and eigenVector_ array are
      //       treated as 2-dimensional arrays of dimension
      //       (numFreq_,2) and (size_,numArnoldiVec_), respectively.
      //       However, taking account of FORTRAN storage scheme for arrays 
      //       (column first) in data transfer, we use a simple 1d vector 
      //       memory allocation here!

      bool rvec = TRUE;
      Double *select = New Double [numArnoldiVec_];
      Double *d = New Double [numArnoldiVec_*2];
      Double omgShift = pow(freqShift_*8.0*atan(1.0),2);

      ARPACK_DSEUPD(&rvec, "A", select, d, eigenVectors_, (Integer*) &size_, 
            &omgShift, "G", (Integer*) &size_, which_, (Integer*) &numFreq_,
            &tolerance_, residual, (Integer*) &numArnoldiVec_, matrixV, 
            (Integer*) &size_, iparams, ipntr, workD, workL, &lenWorkL, &info);

      Integer numEVConverged = iparams[4];
      for (Integer i=0; i< numEVConverged; i++) {
          eigenValues_[i] = d[i];                 // cf. note above!
          // calculate error norm of obtained solution
          vecX = matrixV + i*size_;
          vecY = workD + size_;
          interface_->MultBV(vecX, vecY); 
          vecY = workD;
          interface_->MultAV(vecX, vecY); 
          Double vNrm2 = 0.0;
          for (Integer j=0; j<size_; j++) {
              vecX[j] = vecY[j] - d[i]*vecY[j+size_];
              vNrm2 += vecX[j]*vecX[j];
          }
          eigenTolerances_[i] = sqrt(vNrm2)/d[i];    // cf. note above!
          (*cla) << " Mode no. " << i << ", tol = " << eigenTolerances_[i] << "\n";
      }
    
      return numEVConverged;
  }

  void ArpackSolver::InitTempSpace(Double *tempV,Double *residual,
             Double* workD, Double* workL, Double* matrixV) {


      Integer i;
      for (i=0; i < size_;i++) {
           tempV[i] = 0.0;
           residual[i] = 0.0;
      }

      for (i=0; i < 3*size_;i++)
           workD[i] = 0.0;

      Integer len = numArnoldiVec_*(numArnoldiVec_+8);
      for (i=0; i < len; i++)
          workL[i] = 0.0;

      len = size_*numArnoldiVec_;
      for (i=0; i < len; i++)
           matrixV[i] = 0.0;

  }

  Double ArpackSolver::Eigenvalue(UInt i) {
      return eigenValues_[i];
  }

  Double ArpackSolver::Tolerance(UInt i) {
      return eigenTolerances_[i];
  }

  Double *ArpackSolver::GetEigenvector( UInt modeNr ) {
      return eigenVectors_ + modeNr*size_;
  }

  void ArpackSolver::DebugOff( ) {

      ARPACK_DEBUG_CMN.logfil = 6;
      ARPACK_DEBUG_CMN.ndigit = 0;
      ARPACK_DEBUG_CMN.mgetv0 = 0;
      ARPACK_DEBUG_CMN.msaupd = 0;
      ARPACK_DEBUG_CMN.msaup2 = 0;
      ARPACK_DEBUG_CMN.msaitr = 0;
      ARPACK_DEBUG_CMN.mseigt = 0;
      ARPACK_DEBUG_CMN.msapps = 0;
      ARPACK_DEBUG_CMN.msgets = 0;
      ARPACK_DEBUG_CMN.mseupd = 0;
      ARPACK_DEBUG_CMN.mnaupd = 0;
      ARPACK_DEBUG_CMN.mnaup2 = 0;
      ARPACK_DEBUG_CMN.mnaitr = 0;
      ARPACK_DEBUG_CMN.mneigt = 0;
      ARPACK_DEBUG_CMN.mnapps = 0;
      ARPACK_DEBUG_CMN.mngets = 0;
      ARPACK_DEBUG_CMN.mneupd = 0;
      ARPACK_DEBUG_CMN.mcaupd = 0;
      ARPACK_DEBUG_CMN.mcaup2 = 0;
      ARPACK_DEBUG_CMN.mcaitr = 0;
      ARPACK_DEBUG_CMN.mceigt = 0;
      ARPACK_DEBUG_CMN.mcapps = 0;
      ARPACK_DEBUG_CMN.mcgets = 0;
      ARPACK_DEBUG_CMN.mceupd = 0;

  }

  void ArpackSolver::DebugOn( ) {

      ARPACK_DEBUG_CMN.logfil = 6;
      ARPACK_DEBUG_CMN.ndigit = -3;
      ARPACK_DEBUG_CMN.mgetv0 = 1;
      ARPACK_DEBUG_CMN.msaupd = 2;
      ARPACK_DEBUG_CMN.msaup2 = 2;
      ARPACK_DEBUG_CMN.msaitr = 1;
      ARPACK_DEBUG_CMN.mseigt = 1;
      ARPACK_DEBUG_CMN.msapps = 1;
      ARPACK_DEBUG_CMN.msgets = 1;
      ARPACK_DEBUG_CMN.mseupd = 2;
      ARPACK_DEBUG_CMN.mnaupd = 2;
      ARPACK_DEBUG_CMN.mnaup2 = 2;
      ARPACK_DEBUG_CMN.mnaitr = 1;
      ARPACK_DEBUG_CMN.mneigt = 1;
      ARPACK_DEBUG_CMN.mnapps = 1;
      ARPACK_DEBUG_CMN.mngets = 1;
      ARPACK_DEBUG_CMN.mneupd = 1;
      ARPACK_DEBUG_CMN.mcaupd = 2;
      ARPACK_DEBUG_CMN.mcaup2 = 2;
      ARPACK_DEBUG_CMN.mcaitr = 1;
      ARPACK_DEBUG_CMN.mceigt = 1;
      ARPACK_DEBUG_CMN.mcapps = 1;
      ARPACK_DEBUG_CMN.mcgets = 1;
      ARPACK_DEBUG_CMN.mceupd = 2;

  }

}
