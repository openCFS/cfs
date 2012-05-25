#include <limits>

#include "General/Environment.hh"
#include "General/Exception.hh"

#include "ArpackSolver.hh"

namespace CoupledField {

  ArpackSolver::ArpackSolver() :
    eigenValues_(NULL),
    zEigenValues_(NULL),
    eigenTolerances_(NULL),
    eigenVectors_(NULL),
    zEigenVectors_(NULL)
  {

      shiftAndInvert_ = true;
      freqShift_ = 0.0;
      tolerance_ = 1.e-8;
      maxIterations_ = 5000;
      numArnoldiVec_ = 0;
      interface_ = NULL;
			// cast to char* or receive compiler warning
      which_ = (char*) "LM";
      type_ = (char*) "G";

  }
  
  ArpackSolver::~ArpackSolver() {
    
      delete[] eigenValues_;
      delete[] zEigenValues_;
      delete[] eigenTolerances_;
      delete[] eigenVectors_;
      delete[] zEigenVectors_;
      interface_ = NULL;

  }

  void ArpackSolver::Setup( ArpackMatInterface *matInterface, UInt size, 
                            UInt numFreq, Double freqShift, char* which, 
                            char* type, bool shiftMode ) {

      size_ = size;

      // new reference to matrix interface
      interface_ = matInterface;
    
      // Save frequency parameters
      numFreq_ = numFreq;
      freqShift_ = freqShift;
      shiftAndInvert_ = shiftMode;
      which_ = which;
      type_ = type;

      // set default value for Arnoldi vectors
      numArnoldiVec_ = numFreq_*2;
      
      // check, if number of arnoldi vectors is larger than
      // size of system
      if( numArnoldiVec_ > size_ ) {
        UInt newNumFreq = (UInt) size / 2.0;
        numFreq_ = newNumFreq;
        numArnoldiVec_ = numFreq_*2;
        std::stringstream out;
        WARN( "Number of eigenfrequencies will be re-set to "
              << newNumFreq << " as the number of eigenfrequencies may only be"
              << " 1/2 the number of unknowns of the system");
      }
      

      // eigenvalues, tolerances and vectors
      eigenValues_ = new Double [numArnoldiVec_];
      eigenTolerances_ = new Double [numArnoldiVec_];
      eigenVectors_ = new Double [numFreq_*size_];

  }


  void ArpackSolver::QuadSetup( ArpackMatInterface *matInterface, UInt size, 
              UInt numFreq, Double freqShift, char* which, bool shiftMode ) {

      size_ = size;

      // new reference to matrix interface
      interface_ = matInterface;
    
      // Save frequency parameters
      numFreq_   = numFreq;
      freqShift_ = freqShift;

      shiftAndInvert_ = shiftMode;
      which_ = which;

      // adjust defaults for number of iterations and tolerance
      tolerance_     = 1e-10;
      maxIterations_ = 10000;
      // adjust to higher number of wanted frequencies
      Double logNF = std::log10(numFreq_*1.0);
      if (logNF > 1.0) {
        Integer iterInc = 500*std::floor(5*pow(10,logNF-1));
        maxIterations_ += iterInc;
      }

      // set default value for Arnoldi vectors
      numArnoldiVec_ = numFreq_*2;
      
      // check, if number of arnoldi vectors is larger than
      // size of system
      if( numArnoldiVec_ > size_ ) {
        UInt newNumFreq = (UInt) size / 2.0;
        numFreq_ = newNumFreq;
        numArnoldiVec_ = numFreq_*2;
        std::stringstream out;
        WARN( "Number of eigenfrequencies (" << numFreq << ") has been reset to "
              << newNumFreq << " as the number of eigenfrequencies may only be"
              << " 1/2 the number of unknowns of the system" );
      }
      

      // eigenvalues, tolerances and vectors
      zEigenValues_ = new Complex [numArnoldiVec_];
      eigenTolerances_ = new Double [numArnoldiVec_];
      zEigenVectors_ = new Complex [numFreq_*size_];

      DebugOff();

  }

  void ArpackSolver::SetTolerance(Double tol) {
      //! set ARPACK convergence tolerance in case of user specified value
      tolerance_ = tol;
  }

  void ArpackSolver::SetIterations(Integer maxIt) {
      //! set ARPACK maximum number of iterations in case of user specified value
      maxIterations_ = maxIt;
  }

  void ArpackSolver::SetNumVectors(Integer numVec) {
      //! set ARPACK number of Arnoldi vectors in case of user specified value
      numArnoldiVec_ = numVec;
  }

  UInt ArpackSolver::GetNev( ) {
      return numFreq_;
  }

  char* ArpackSolver::GetWhich( ) {
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

      bool converged = false;
      Integer  ido = 0, info = 0;

      // temp vector to store B*x
      Double *tempV = new Double [size_];
      // pointers to the position of x and y in workD
      Double *vecX, *vecY;

      // temp working space required in dsaupd
      Double *workD = new Double [3*size_];
      Double *residual = new Double [size_];
      Double *workL = new Double [numArnoldiVec_*(numArnoldiVec_+8)];
      Double *matrixV = new Double [size_*numArnoldiVec_];

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
      for (itNum=0; itNum<maxIterations_; itNum++) {

          // cast to char* or receive compiler warning
          F77NAME(dsaupd)(&ido, type_, (Integer*) &size_, which_, (Integer*) &numFreq_, 
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
            EXCEPTION( "User required shifts not yet implemented!\n"
                       << "ARPACK requested ido=3 in DSAUPD!" );
            break;

          case 99:
            converged = true;
            break;

          default:
              // something went wrong -> arpack requires undefind ido
            EXCEPTION( "ARPACK requested not implemented ido=" << ido 
                       << " in DSAUPD!" );
            break;

          }

          // if convergence has been achieved, exit loop
          if (converged)
              break;

      }

      // check whether everything went well
      if (info<0) {
          EXCEPTION("Error reported in ritz value calculation:\n"
           << ArpackError(info) );
      }

      if ( itNum>maxIterations_ && info != 99 ) {
        EXCEPTION("Maximum number of iterations has been reached, but no solution"
                 << "could be obtained which satisfies requested tolerance setting.\n"
                 << "Increase maxIt (=" << maxIterations_ << ") in xml file!" );
      }

      // proceed with calculation of eigenvectors, values, and tolerances

      // NOTE: in the FORTRAN routine dseupd d and eigenVector_ array are
      //       treated as 2-dimensional arrays of dimension
      //       (numFreq_,2) and (size_,numArnoldiVec_), respectively.
      //       However, taking account of FORTRAN storage scheme for arrays 
      //       (column first) in data transfer, we use a simple 1d vector 
      //       memory allocation here!

      bool rvec = true;
      Double *select = new Double [numArnoldiVec_];
      Double *d = new Double [numArnoldiVec_*2];
      Double omgShift = pow(freqShift_*8.0*atan(1.0),2);

			// cast to char* or receive compiler warning
      F77NAME(dseupd)(&rvec, (char*) "A", select, d, eigenVectors_, (Integer*) &size_, 
                      &omgShift, (char*) "G", (Integer*) &size_, which_, (Integer*) &numFreq_,
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
          for (UInt j=0; j<size_; j++) {
              vecX[j] = vecY[j] - d[i]*vecY[j+size_];
              vNrm2 += vecX[j]*vecX[j];
          }
          eigenTolerances_[i] = sqrt(vNrm2)/d[i];    // cf. note above!
          //(*cla) << " Mode no. " << i << ", tol = " << eigenTolerances_[i] << "\n";
      }
    
      return numEVConverged;
  }

  UInt ArpackSolver::FindQuadEigenvalues( ) {

      bool converged = false;
      Integer  ido = 0, info = 0;

      // temp vector to store B*x
      Complex *tempV = new Complex [size_];
      // pointers to the position of x and y in workD
      Complex *vecX, *vecY;

      // temp working space required in znaupd
      Complex *workD = new Complex [3*size_];
      Complex *residual = new Complex [size_];
      Integer lenWorkL = numArnoldiVec_*(3*numArnoldiVec_+5);
      Complex *workL = new Complex [lenWorkL];
      Complex *matrixV = new Complex [size_*numArnoldiVec_];
      Double  *workDbleD = new Double [numArnoldiVec_];

      InitQuadTempSpace(tempV, residual, workD, workL, matrixV, workDbleD);

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

        F77NAME(znaupd)(&ido, (char*) "G", (Integer*) &size_, which_,
                        (Integer*) &numFreq_, &tolerance_, residual,
                        (Integer*) &numArnoldiVec_, matrixV,
                        (Integer*) &size_, iparams, ipntr, workD, workL,
                        &lenWorkL, workDbleD, &info);

          switch (ido)  {

          case -1:
              // NOTE: addresses in workd are given as FORTRAN addreses
              // solve (A+shift*D+shift**2*B)*y = B*x
              // x is returned by znaupd in workd(ipntr[1]:ipntr[1]+size-1)
              // y is expected by znaupd in workd(ipntr[2]:ipntr[2]+size-1)

             vecX = workD + --ipntr[0];
             vecY = workD + --ipntr[1];
             interface_->MultBV(vecX, tempV); 
             interface_->MultShiftOpV(tempV, vecY); 
             break;

          case 1:
              // solve (A+shift*D+shift**2*B)*y = B*x
              // B*x is returned by znaupd in workd(ipntr[3]:ipntr[3]+size-1)
              // y is expected by znaupd in workd(ipntr[2]:ipntr(2]+size-1)
             vecX = workD + --ipntr[2];
             vecY = workD + --ipntr[1];
             interface_->MultShiftOpV(vecX,vecY);
             break;

          case 2:
              // calculate y = B*x
              // x is returned by znaupd in workd(ipntr[1]:ipntr[1]+size-1)
              // y is expected by znaupd in workd(ipntr[2]:ipntr[2]+size-1)
              vecX = workD + --ipntr[0];
              vecY = workD + --ipntr[1];
              interface_->MultBV(vecX,vecY);
              break;

          case 3:
              // algorithm requires shifts - > not implemented
              EXCEPTION("User required shifts not yet implemented!\n"
                       << "ARPACK requested ido=3 in ZNAUPD!");
              break;

          case 99:
              converged = true;
              break;

          default:
              // something went wrong -> arpack requires undefind ido
              EXCEPTION("ARPACK requested not implemented ido=" << ido 
                       << " in ZNAUPD!");
              break;

          }

          // if convergence has been achieved, exit loop
          if (converged)
              break;


      }

      // check whether everything went well
      if (info != 0) {
        EXCEPTION("iparam[5] " << iparams[4] << 
        "Error reported in ritz value calculation:\n" << ArpackError(info) );
      }

      if ( itNum>maxIterations_ && info != 99 ) {
        EXCEPTION("Maximum number of iterations has been reached, but no solution"
            << "could be obtained which satisfies requested tolerance setting.\n"
            << "Increase maxIt (=" << maxIterations_ << ") in xml file!");
      }

      // proceed with calculation of eigenvectors, values, and tolerances

      // NOTE: in the FORTRAN routine zneupd d and eigenVector_ array are
      //       treated as 2-dimensional arrays of dimension
      //       (numFreq_,2) and (size_,numArnoldiVec_), respectively.
      //       However, taking account of FORTRAN storage scheme for arrays 
      //       (column first) in data transfer, we use a simple 1d vector 
      //       memory allocation here!

      bool rvec = true;
      Double *select = new Double [numArnoldiVec_];
      Complex *d = new Complex [numArnoldiVec_*2];
      Complex *zwv = new Complex [numArnoldiVec_*2];

      // we define a second tmp-array to store all eigenvectors, 
      // also for the negative frequencies. 
      // zEigenVectors_ will only store data for positive freqs.
      Complex *zEV = new Complex [numFreq_*size_];

      for (UInt i=0; i< numArnoldiVec_*2; i++) {
        d[i] = (Complex) 0.0;
	zwv[i] = (Complex) 0.0;
      }
      for (UInt i=0; i< numFreq_*size_; i++) {
        zEV[i]            = (Complex) 0.0;
        zEigenVectors_[i] = (Complex) 0.0;
      }

      Complex omgShift = (Complex) freqShift_*8.0*atan(1.0);

      F77NAME(zneupd)(&rvec, (char*) "A", select, d, zEV,
                      (Integer*) &size_, &omgShift, zwv, (char*) "G",
                      (Integer*) &size_, which_, (Integer*) &numFreq_,
                      &tolerance_, residual, (Integer*) &numArnoldiVec_,
                      matrixV, (Integer*) &size_, iparams, ipntr, workD,
                      workL, &lenWorkL, workDbleD, &info);

      Integer numEVConverged = iparams[4];

      Complex *z1 = new Complex [size_/2];
      Complex *z2 = new Complex [size_/2];
      Complex *z3 = new Complex [size_/2];
      Complex *zA = new Complex [size_];
      Complex *zB = new Complex [size_];
      Complex *zC = new Complex [size_];
      Complex *ztmp;

      // we only count and treat/store eigenvalues with positive frequency value
      Integer *sortedEV = new Integer [numEVConverged];
      Double  *sortedFreq = new Double [numEVConverged];
      for (Integer ncv=0; ncv < numEVConverged; ncv++) {
        sortedEV[ncv] = -1;
        sortedFreq[ncv] = 0.0;
      }
      
      Integer found = 0, pos;
      for (Integer ncv=0; ncv < numEVConverged; ncv++) {
        Double freq = d[ncv].imag();
        if (freq >= 0.0) {
          if (found==0) {
             sortedFreq[found] = freq;
             sortedEV[found] = ncv;
          } else {
            pos = found;
            for (Integer k=found-1; k >=0; k--) {
              if (sortedFreq[k] > freq)
                 pos--;
            }
            for (Integer k=found-1; k >=pos; k--) {
              sortedFreq[k+1] = sortedFreq[k];
              sortedEV[k+1] = sortedEV[k];
            }
            sortedFreq[pos] = freq;
            sortedEV[pos] = ncv;
          }
          found++;
        }
      }
//      std::cout << "sortedEV, sortedFreq " << std::endl;
//      for (Integer pos=0; pos<numEVConverged; pos++) {
//        std::cout << "pos " << pos << " nev " << sortedEV[pos] << ", freq " 
//                  << sortedFreq[pos] << std::endl;
//      }

      Integer *evPos = new Integer [numEVConverged];
      for (Integer pos=0; pos<numEVConverged; pos++) {
        evPos[pos] = -1;
      }
      for (Integer pos=0; pos<found; pos++) {
	if (sortedEV[pos] >= 0) {
          evPos[sortedEV[pos]] = pos;
        }
      }
      
//      std::cout << "evPos " << std::endl;
//      for (Integer pos=0; pos<numEVConverged; pos++) {
//        std::cout << "pos " << pos << " nev " << evPos[pos] << std::endl;
//      }

      for (Integer ncv=0; ncv < numEVConverged; ncv++) {

        if (d[ncv].imag() >= 0.0) {

          pos = evPos[ncv];
	  if (pos >= 0) {
            zEigenValues_[pos] = d[ncv];
            for (UInt j=0; j<size_; j++) {
              zEigenVectors_[pos*size_+j] = zEV[ncv*size_+j];
            }

            ztmp = zEigenVectors_ + pos*size_;
            for (UInt j=0; j<size_/2; j++) {
              z3[j] = ztmp[j];
            }
            interface_->MultZDampV(ztmp, z1); 
            ztmp += size_/2;
            interface_->MultZStiffV(ztmp, z2); 
            for (UInt j=0; j<size_/2; j++) {
              z1[j] = -z1[j] - z2[j];
            }
  
            if ( interface_->UseStiffInNMat() ) {
              for (UInt j=0; j<size_/2; j++) {
                z2[j] = z3[j];
              }
              interface_->MultZStiffV(z2, z3); 
            }
            interface_->ScaleDiag(z3);
            for (UInt j=0; j<size_/2; j++) {
                zA[j]         = z1[j];
                zA[j+size_/2] = z3[j];
            }
  
            ztmp = zEigenVectors_ + pos*size_;
            for (UInt j=0; j<size_/2; j++) {
               z3[j] = ztmp[j+size_/2];
            }
            interface_->MultZMassV(ztmp, z1); 
  
            if ( interface_->UseStiffInNMat() ) {
              for (UInt j=0; j<size_/2; j++) {
                  z2[j] = z3[j];
              }
              interface_->MultZStiffV(z2, z3); 
            }
            interface_->ScaleDiag(z3);
            for (UInt j=0; j<size_/2; j++) {
              zB[j]         = z1[j];
              zB[j+size_/2] = z3[j];
            }
  
            Double vNrm2 = 0.0,v2Nrm2 = 0.0;
            for (UInt j=0; j<size_/2; j++) {
              zC[j] = zA[j]-d[ncv]*zB[j];
              vNrm2 += pow(abs(zC[j]),2);
            }
            for (UInt j=size_/2; j<size_; j++) {
              zC[j] = zA[j]-d[ncv]*zB[j];
              v2Nrm2 += pow(abs(zC[j]),2);
            }
            vNrm2 += v2Nrm2;
            eigenTolerances_[pos] = sqrt(vNrm2)/abs(d[ncv]);
  
            // "normalize" eigenvectors: 
            ztmp = zEigenVectors_ + pos*size_;
            Double  epsTest = 1.0e-6;
            // 1) for a unit vector
            Complex vLen = (Complex) 0.0;
            for (UInt j=size_/2; j<size_; j++) {
              vLen += ztmp[j]*std::conj(ztmp[j]);
            }
            vLen = sqrt(vLen);
            for (UInt j=0; j<size_; j++) {
              ztmp[j] /= vLen;
            }
  
            // 2) start at size_/2, i.e. use "essential part" of the eigenvector
            //    and look for first non-noisy entry
            UInt jNorm = size_/2, jFound=0;
            while (jFound==0 && jNorm<size_) {
              if (std::abs(ztmp[jNorm]) > epsTest)
                jFound = 1;
              else
                jNorm++;
            }
  
            // 3) rescale complete vector so that found entry will be scaled
            //    with equal real and imag parts. this is done to get rid of
            //    numerical noise in comparing eigenvectors
            //    NOTE: this should be modified, as soon as complex valued 
            //    eigenvectors are available for the result files
  
            Complex zbase = Complex (1.,1.);
            Complex zNorm = ztmp[jNorm]*zbase;
            for (UInt j=0; j<size_; j++) {
              ztmp[j] /= zNorm;
            }

	    // 4) calculated conj(z)^T * M * z
            //          ztmp = zEigenVectors_ + pos*size_ + size_/2;
            //          interface_->MultZMassV(ztmp, z1); 
            //          Complex normalize = (Complex) 0.0;
            //          for (Integer j=0; j<size_/2; j++) {
            //             normalize += std::conj(ztmp[j])*z1[j];
            //          }
            //	  normalize = sqrt(normalize);
            //          ztmp = zEigenVectors_ + pos*size_;
            //          for (Integer j=0; j<size_; j++) {
            //             ztmp[j] /= normalize;
            //          }
	    //-------------------------------------
          }
        }

      }
      for (pos=0; pos<found; pos++) {
          (*cla) << " EV no. " << pos << ", EV  = " << zEigenValues_[pos] << "\n";
          (*cla) << " EV no. " << pos << ", tol = " << eigenTolerances_[pos] << "\n";
      }
    
      return found;
  }

  void ArpackSolver::InitTempSpace(Double *tempV,Double *residual,
             Double* workD, Double* workL, Double* matrixV) {


      Integer i;
      for (i=0; i < (int) size_;i++) {
           tempV[i] = 0.0;
           residual[i] = 0.0;
      }

      for (i=0; i < 3*((int) size_);i++)
           workD[i] = 0.0;

      Integer len = numArnoldiVec_*(numArnoldiVec_+8);
      for (i=0; i < len; i++)
          workL[i] = 0.0;

      len = size_*numArnoldiVec_;
      for (i=0; i < len; i++)
           matrixV[i] = 0.0;

  }

  void ArpackSolver::InitQuadTempSpace(Complex *tempV,Complex *residual,
         Complex* workD, Complex* workL, Complex* matrixV, Double* workDbleD) {

      UInt i;
      for (i=0; i < size_;i++) {
           tempV[i] = (Complex) 0.0;
           residual[i] = (Complex) 0.0;
      }

      for (i=0; i < 3*size_;i++)
           workD[i] = (Complex) 0.0;

      UInt len = numArnoldiVec_*(3*numArnoldiVec_+5);
      for (i=0; i < len; i++)
          workL[i] = (Complex) 0.0;

      len = size_*numArnoldiVec_;
      for (i=0; i < len; i++)
           matrixV[i] = (Complex) 0.0;

      for (i=0; i < numArnoldiVec_; i++)
           workDbleD[i] = 0.0;
  }

  Double ArpackSolver::Eigenvalue(UInt i) {
      return eigenValues_[i];
  }

  Complex ArpackSolver::CmplxEigenvalue(UInt i) {
      return zEigenValues_[i];
  }

  Double ArpackSolver::Tolerance(UInt i) {
      return eigenTolerances_[i];
  }

  Double *ArpackSolver::GetEigenvector( UInt modeNr ) {
      return eigenVectors_ + modeNr*size_;
  }
  
  Complex *ArpackSolver::GetQuadEigenvector( UInt modeNr ) {
      return zEigenVectors_ + modeNr*size_;
  }
  
  std::string ArpackSolver::ArpackError( Integer errNo ) {
    std::string msg;

    switch (errNo ) {
      case 0:
        msg = "No error";
        break;

      case 1: 
        msg = "Maximum number of iterations taken. "
          " All possible eigenvalues of OP has been found. IPARAM(5) "
          " returns the number of wanted converged Ritz values";
        break;

      case 2: 
        msg = "No longer an informational error. Deprecated starting"
          "with release 2 of ARPACK.";
        break;

      case 3: 
        msg = "No shifts could be applied during a cycle of the"
          "Implicitly restarted Arnoldi iteration. One possibility"
          "is to increase the size of NCV relative to NEV";
        break;

      case  -1: 
        msg = "N must be positive";
        break;

      case -2: 
        msg = "NEV must be positive";
        break;

      case -3: 
        msg = "NCV must be greater than NEV and less than or equal to N";
        break;
      
      case -4: 
        msg = "The maximum number of Arnoldi update iterations allowed"
        "must be greater than zero.";
        break;
      
      case -5: 
        msg = "WHICH must be one of 'LM', 'SM', 'LA', 'SA' or 'BE'";
        break;
      
      case -6: 
        msg = "BMAT must be one of 'I' or 'G'";
        break;
      
      case -7: 
        msg = "Length of private work array WORKL is not sufficient.";
        break;
      
      case -8: 
        msg = "Error return from trid. eigenvalue calculation;" 
          " Informatinal error from LAPACK routine dsteqr.";
        break;
      
      case -9: 
        msg = "Starting vector is zero.";
        break;
      
      case -10: 
        msg = "IPARAM(7) must be 1,2,3,4,5";
        break;
      
      case -11: 
        msg = "IPARAM(7) = 1 and BMAT = 'G' are incompatable.";
        break;
      
      case -12: 
        msg = "IPARAM(1) must be equal to 0 or 1.";
        break;
      
      case -13: 
        msg = "NEV and WHICH = 'BE' are incompatable.";
        break;
      
      case -9999: 
      msg = "Could not build an Arnoldi factorization."
        " IPARAM(5) returns the size of the current Arnoldi"
        " factorization. The user is advised to check that"
        " enough workspace and array storage has been allocated.";
        
        break;

      default:
        msg = "An unknown error occured";
    }
    return msg;

  }

  void ArpackSolver::DebugOff( ) {

      ::F77NAME(debug).logfil = 6;
      ::F77NAME(debug).ndigit = 0;
      ::F77NAME(debug).mgetv0 = 0;
      ::F77NAME(debug).msaupd = 0;
      ::F77NAME(debug).msaup2 = 0;
      ::F77NAME(debug).msaitr = 0;
      ::F77NAME(debug).mseigt = 0;
      ::F77NAME(debug).msapps = 0;
      ::F77NAME(debug).msgets = 0;
      ::F77NAME(debug).mseupd = 0;
      ::F77NAME(debug).mnaupd = 0;
      ::F77NAME(debug).mnaup2 = 0;
      ::F77NAME(debug).mnaitr = 0;
      ::F77NAME(debug).mneigt = 0;
      ::F77NAME(debug).mnapps = 0;
      ::F77NAME(debug).mngets = 0;
      ::F77NAME(debug).mneupd = 0;
      ::F77NAME(debug).mcaupd = 0;
      ::F77NAME(debug).mcaup2 = 0;
      ::F77NAME(debug).mcaitr = 0;
      ::F77NAME(debug).mceigt = 0;
      ::F77NAME(debug).mcapps = 0;
      ::F77NAME(debug).mcgets = 0;
      ::F77NAME(debug).mceupd = 0;

  }

  void ArpackSolver::DebugOn( ) {

      ::F77NAME(debug).logfil = 6;
      ::F77NAME(debug).ndigit = -3;
      ::F77NAME(debug).mgetv0 = 1;
      ::F77NAME(debug).msaupd = 2;
      ::F77NAME(debug).msaup2 = 2;
      ::F77NAME(debug).msaitr = 1;
      ::F77NAME(debug).mseigt = 1;
      ::F77NAME(debug).msapps = 1;
      ::F77NAME(debug).msgets = 1;
      ::F77NAME(debug).mseupd = 2;
      ::F77NAME(debug).mnaupd = 2;
      ::F77NAME(debug).mnaup2 = 2;
      ::F77NAME(debug).mnaitr = 1;
      ::F77NAME(debug).mneigt = 1;
      ::F77NAME(debug).mnapps = 1;
      ::F77NAME(debug).mngets = 1;
      ::F77NAME(debug).mneupd = 1;
      ::F77NAME(debug).mcaupd = 2;
      ::F77NAME(debug).mcaup2 = 2;
      ::F77NAME(debug).mcaitr = 1;
      ::F77NAME(debug).mceigt = 1;
      ::F77NAME(debug).mcapps = 1;
      ::F77NAME(debug).mcgets = 1;
      ::F77NAME(debug).mceupd = 2;

  }

}
